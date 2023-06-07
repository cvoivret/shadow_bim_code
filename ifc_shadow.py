from collections import defaultdict
import itertools
from array import array
import numpy as np
import pandas as pd
import ifcopenshell

from timeit import default_timer as timer

from ifcopenshell.geom import create_shape
from ifcopenshell.geom.occ_utils import yield_subshapes

from OCC.Display.SimpleGui import init_display
from OCC.Core.Quantity import Quantity_Color,Quantity_TOC_RGB

from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox,BRepPrimAPI_MakePrism,BRepPrimAPI_MakeHalfSpace,BRepPrimAPI_MakeSphere,BRepPrimAPI_MakeCylinder
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties,brepgprop_VolumeProperties,brepgprop_LinearProperties
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Sewing,BRepBuilderAPI_MakeSolid
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepTools import breptools_UVBounds

from OCC.Core.GProp import GProp_GProps
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.gp import gp_Pnt,gp_Dir,gp_Vec,gp_Pln,gp_Lin,gp_Trsf,gp_Ax3
from OCC.Core.Geom import Geom_Plane

from OCC.Core.TopoDS import TopoDS_Face
from OCC.Core.TopTools import TopTools_ListOfShape,TopTools_IndexedMapOfShape
from OCC.Core.TopExp import topexp_MapShapes
from OCC.Core.TopAbs import TopAbs_SOLID,TopAbs_FACE,TopAbs_SHELL,TopAbs_WIRE

from OCC.Extend.TopologyUtils import TopologyExplorer, WireExplorer
 
from OCC.Core.BOPAlgo import BOPAlgo_BOP,BOPAlgo_Operation
from OCC.Core.BOPAlgo import BOPAlgo_CellsBuilder
from OCC.Core.BOPTools import BOPTools_AlgoTools_OrientFacesOnShell

from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib

import OCC.Core.ShapeFix as ShapeFix_Shape

from OCC.Core.IntCurvesFace import IntCurvesFace_ShapeIntersector

from OCC.Extend.DataExchange import write_stl_file


def fuse_listOfShape(los,FuzzyValue=1e-6):
    """
    Simple funtion to wrap boilerplate code
    """    
    fuser=BOPAlgo_BOP()
    fuser.SetOperation(BOPAlgo_Operation.BOPAlgo_FUSE)
    los_1=TopTools_ListOfShape()
    [los_1.Append(s) for s in los[::2]]
    los_2=TopTools_ListOfShape()
    [los_2.Append(s) for s in los[1::2]]
    fuser.SetArguments(los_1)
    fuser.SetTools(los_2)
    fuser.SetFuzzyValue(FuzzyValue)
    fuser.SetNonDestructive(True)
    fuser.Perform()
    return fuser.Shape()

def shapes_as_solids(lshape):
    """
    Try to build a list of solid from a list to shapes.
    Flatten nested shapes.
    Sew tesselated shape to build a solid.
    Fix shape if needed.
    
    """
    lsolid=[]
     
    maps=TopTools_IndexedMapOfShape()
    for s in lshape:
        maps.Clear()
        topexp_MapShapes(s,TopAbs_SOLID,maps)
        if(maps.Size()>0):
            lsolid.extend([maps.FindKey(i) for i in range(1,maps.Size()+1)])
        else:
            maps.Clear()
            topexp_MapShapes(s,TopAbs_FACE,maps)
            sewer=BRepBuilderAPI_Sewing()
            [sewer.Add(maps.FindKey(i)) for i in range(1,maps.Size()+1)]
            sewer.Perform()
            sewed=sewer.SewedShape()
            if(sewed.ShapeType()==0):
                lshell=list(yield_subshapes(sewed))
                
                for shell in lshell:
                    lsolid.append(BRepBuilderAPI_MakeSolid(shell).Solid())
            else:
                solid=BRepBuilderAPI_MakeSolid(sewed).Solid()
                lsolid.append(solid)
    lsolid2=[]            
    for s in lsolid:
        fixer=ShapeFix_Shape.ShapeFix_Shape(s)
        fixer.Perform()
        lsolid2.append(fixer.Shape())
         
    return lsolid2

def get_external_shell(lshape):
    """
    try to identigy a shell (set of face) that limit the inside and the outside of the building.
    Basically, wall and room must be part of the solid list input.
    Build the boundign box of the whole model and enlarge it.
    
    """
       
    #unionize solids
    unionsolid=fuse_listOfShape(lshape)
    
    # Create the bounding box of the unioned model
    box=Bnd_Box()
    brepbndlib.Add(unionsolid,box)
    box.Enlarge(1.) # to avoid parallel face of bbox with face model to be coplanar(
    # broke shell intersection 
    boxshape=BRepPrimAPI_MakeBox(box.CornerMin(),box.CornerMax()).Shape()

    #boolean difference between the unioned model and its bounding box
    diff=BOPAlgo_BOP()
    diff.SetOperation(BOPAlgo_Operation.BOPAlgo_CUT)
    diff.AddArgument(boxshape)
    diff.AddTool(unionsolid)
    diff.SetFuzzyValue(1e-5)
    diff.Perform()
    diffshape=diff.Shape()
    

    # boolean common of shells : could be considered as the shell 
    # separating interior and exterior of the building
    top=TopologyExplorer(unionsolid)
    unionshell = top.shells()
    

    top=TopologyExplorer(diffshape)
    diffshell = top.shells()


    common=BOPAlgo_BOP()
    common.SetOperation(BOPAlgo_Operation.BOPAlgo_COMMON)
    args=TopTools_ListOfShape()
    [args.Append(shell) for shell in unionshell]
    
    tools=TopTools_ListOfShape()
    [tools.Append(shell) for shell in diffshell]
    common.SetArguments(args)
    common.SetTools(tools)
    common.SetFuzzyValue(1e-5)
    common.Perform()
    commonshell=common.Shape()
    
    BOPTools_AlgoTools_OrientFacesOnShell(commonshell)
    
    return commonshell
 

def shadow_caster_ext(sun_dir,building,theface,theface_norm,min_area = 1e-3):
    """
    sun_dir = one vector (downward direction)
    building = a solids that possibily make shadow on face
    face = a face to cast shadow on from building along sun_dir
    face_norm = pointing to the exterior of the face (outside)
    
    return  : a face with zero or positive area, None if no shadow
    
    """
    #print(theface_norm.Dot(sun_dir))
    # face not exposed to the sun
    if theface_norm.Dot(sun_dir)>-1e-5:
        #print('not exposed',flush=True)
        return theface# void face with zero area
    gpp=GProp_GProps()
    brepgprop_SurfaceProperties(theface,gpp)
    gf_area=gpp.Mass()
    
    ext_vec=gp_Vec(sun_dir)
    ext_vec.Multiply(5)
    
    # extrusion of 
    extrusion1=BRepPrimAPI_MakePrism(theface,-ext_vec,False,True).Shape()
    
    intersector=BOPAlgo_BOP()
    intersector.SetOperation(BOPAlgo_Operation.BOPAlgo_COMMON)
    intersector.AddTool(extrusion1) 
    intersector.AddArgument(building)
    intersector.Perform()
    intersection=intersector.Shape()
        
    intersection_faces=list(TopologyExplorer(intersection).faces())
               
    larea=[]
    lfaces=[]
    
    for ff in intersection_faces:
        
        adapt=BRepAdaptor_Surface(ff)
        if adapt.GetType()==1:
            cyl=adapt.Cylinder()
            umin,umax,vmin,vmax=breptools_UVBounds(ff)
            if vmin<0.0:
                cyl.VReverse()
            
            ax3=cyl.Position()
            vec=gp_Dir(*sun_dir.Coord())
            
            vec.Cross(ax3.Direction())
            newax3=gp_Ax3(ax3.Location(),ax3.Direction(),vec)
            shape=BRepPrimAPI_MakeCylinder(newax3.Ax2(),cyl.Radius()*2,2,3.14).Shape()
            
            com=BRepAlgoAPI_Common(shape,ff)
            com.Build()
            shape=com.Shape()
            #lcyl.append(shape)
            maps=TopTools_IndexedMapOfShape()
            topexp_MapShapes(shape,TopAbs_FACE,maps)
            lfacetokeep=[maps.FindKey(i) for i in range(1,maps.Size()+1)]
            if( len(lfacetokeep)==1):
                ff=lfacetokeep[0]
            else:
                continue
        
        srf3 = BRep_Tool().Surface(ff)
        umin,umax,vmin,vmax=breptools_UVBounds(ff)
        props=GeomLProp_SLProps(srf3,0.5*(umax-umin),0.5*(vmax-vmin),1,0.001)
        fn=props.Normal()
        
        
        
        if(ff.Orientation()==1):
            fn.Reverse()
        # avoid face nearly parallel with extrusion generatrix
        # ie face with normal perpendicular with extrusion direction
        if(fn.Dot(sun_dir)<-1e-5):
            brepgprop_SurfaceProperties(ff,gpp)
            larea.append(gpp.Mass())
            if(ff.Orientation()==1):
                ff.Reverse()
            
            lfaces.append(ff)
    
    lsolid=[ BRepPrimAPI_MakePrism(s,ext_vec,False,True).Shape() for s in lfaces]
    
    
    if(len(lsolid)==0):
        return TopoDS_Face() # void face with zero area
    
    brepgprop_SurfaceProperties(theface,gpp)
    totarea=gpp.Mass()
    
    lface2=[]
    for s,f in zip(lsolid,lfaces):
        common=BRepAlgoAPI_Common(s,theface)
        common.Build()
        sh=common.Shape()
        
        brepgprop_SurfaceProperties(sh,gpp)
        area_proj=gpp.Mass()
        #brepgprop_SurfaceProperties(f,gpp)
        #area=gpp.Mass()
        if(area_proj/totarea<1e-4):
            continue
        lface2.append(sh)
        
    
    if len(lface2)==1:
        shadowface=lface2[0]
    else:    
        los2 = TopTools_ListOfShape()
        [los2.Append(s) for s in lface2]
        
        cb=BOPAlgo_CellsBuilder()
        cb.SetArguments(los2)
        cb.Perform()
        cb.AddAllToResult(2,False)
        cb.RemoveInternalBoundaries()
        shadowface=cb.Shape()
    
    if shadowface==None:
        return TopoDS_Face()
       
    return shadowface




def shadow_caster_ray(sun_dir,building,theface,theface_norm,Nray=5):
 
    sphere_rad=0.05
    lshape=[]
        
    #discretize the face with Nray points
    srf = BRep_Tool().Surface(theface)
    umin,umax,vmin,vmax=breptools_UVBounds(theface)
    
    uoffset=0.5*(umax-umin)/Nray
    voffset=0.5*(vmax-vmin)/Nray
    
    uvalues,vvalues= np.meshgrid(np.linspace(umin+uoffset,umax-uoffset,Nray),
                                 np.linspace(vmin+voffset,vmax-voffset,Nray))
    
    # face not exposed to the sun
    if theface_norm.Dot(sun_dir)>-1.e-5:
        
        for u,v in zip(uvalues.flatten(),vvalues.flatten()):
            point=srf.Value(u,v)
            #lshape.append(BRepPrimAPI_MakeSphere(point,sphere_rad).Shape())
        return np.ones(uvalues.shape)#,lshape# all points of discretization are in shadow
    
    
    shape_inter = IntCurvesFace_ShapeIntersector()
    shape_inter.Load(building, 1e-6)
    infinity=float("+inf")
    nbpoints=array('b')
    for u,v in zip(uvalues.flatten(),vvalues.flatten()):
        point=srf.Value(u,v)
        line=gp_Lin(point,-sun_dir)
        shape_inter.PerformNearest(line, 0.0,100.)
        nbpoints.append(shape_inter.NbPnt())
        #if(shape_inter.NbPnt()>0):
        #    lshape.append(BRepPrimAPI_MakeSphere(point,sphere_rad).Shape())
    
    #print(nbpoints)
    res=np.array(nbpoints).reshape(uvalues.shape)
    
    res[res>0.]=1
    return res 

def exterior_wall_normal(wall_shape,external_shell):
    
    gpp=GProp_GProps()
    
    tools=TopTools_ListOfShape()
    tools.Append(external_shell) 

    args=TopTools_ListOfShape()
    args.Append(wall_shape)
        
    common=BOPAlgo_BOP()
    common.SetOperation(BOPAlgo_Operation.BOPAlgo_COMMON)
        
    common.SetArguments(args)
    common.SetTools(tools)
    
    common.SetFuzzyValue(1e-6)
    common.Perform()
    commonshell2=common.Shape() 
        
    
    # exteriro wall !!
    if commonshell2:
        faces=list(TopologyExplorer(commonshell2).faces())
        norm_area=defaultdict(float)
        norm_map=defaultdict(list)
        for f in faces:
            srf = BRep_Tool().Surface(f)
            plane = Geom_Plane.DownCast(srf)
            fn = plane.Axis().Direction()
            if(f.Orientation()==1):
                fn.Reverse()
            face_norm_coord=fn.Coord()
            # maybe necessary to round...
            face_norm_coord = tuple(round(c,10) for c in face_norm_coord)
            brepgprop_SurfaceProperties(f,gpp)
            norm_area[face_norm_coord]+=gpp.Mass()
            norm_map[face_norm_coord].append(f)
        wall_norm = max(norm_area, key=norm_area.get)   
                
        #print(norm_area)
        #print(wall_norm)
        #print(norm_map[wall_norm])
                    
        # wall_norm is rounded but almost equal to all element in the list
        # taking the first
        #lface_wall.append(norm_map[wall_norm])
        #lface_wall.append(norm_map[wall_norm][0])
        first_wall_face =norm_map[wall_norm][0]
        srf = BRep_Tool().Surface(first_wall_face)
        plane = Geom_Plane.DownCast(srf)
        wall_norm = plane.Axis().Direction()
        if(first_wall_face.Orientation()==1):
            wall_norm.Reverse()
        
        return wall_norm
    

def exterior_wall_normal_dict(wallwindow,external_shell):
    
    wallnorm=dict()
    #shape of wall with a window
    wall_shapes=[create_shape(setting, ifc_file.by_guid(w_id)).geometry 
                    for w_id in wallwindow.keys() if ifc_file.by_guid(w_id).Representation 
                    is not None]
    for (w_id,ws) in zip(wallwindow.keys(),wall_shapes):
        wallnorm[w_id]=exterior_wall_normal(ws,external_shell)
    
    return wallnorm
    
def biggestface_along_vector(shape,vector,tol=1e-6,ratio=0.9):
    gpp=GProp_GProps()
    faces=list(TopologyExplorer(shape).faces())
    #print(" nb face par fenetre ", len(faceswin))
    facelist=[]
    facearea=[]
    #facenormal=[]
    for f in faces:
        top=TopologyExplorer(f)
        #print(top.number_of_wires())
        # face with some kind of hole
        if top.number_of_wires()>1:
            continue
        srf = BRep_Tool().Surface(f)
        plane2 = Geom_Plane.DownCast(srf)
        face_norm = plane2.Axis().Direction()
        if(f.Orientation()==1):
            face_norm.Reverse()
        
        
        if(face_norm.IsEqual(vector,tol)):
            #print(" face2 ",win_norm.Coord())
           
            brepgprop_SurfaceProperties(f,gpp)
            #print(" area ", gpp.Mass())
            facearea.append(round(gpp.Mass(),5))
            facelist.append(f)
            #facenormal.append(face_norm)
    #print('\n window ',i)
    
    maxarea=max(facearea)
    gfaces=[ face for area,face in zip(facearea,facelist) if 
                area>maxarea*ratio]
    return gfaces

def biggestfaces_along_normaldict(wallwindow,wallnormal):
    glassface_bywindowid=defaultdict(list)
    #gpp=GProp_GProps()    
    #print(" wall norm ", wall_norm.Coord())
    for w_id in wallwindow.keys():
        if (w_id in wallnorm.keys()):
            wall_norm=wallnormal[w_id]
        else:
            # window in interior wall
            continue
            
        for win_id in wallwindow[w_id]:
        
            windowshape=create_shape(setting, ifc_file.by_guid(win_id)).geometry
            gfaces=biggestface_along_vector(windowshape,wall_norm)
            glassface_bywindowid[win_id].extend(gfaces)
            
    return glassface_bywindowid
 
def window_in_wall(ifcwall):
    windows=[]
    for op in ifcwall.HasOpenings:
            #print('\n ***',w)
            #print('  ',op)
            for re in op.RelatedOpeningElement.HasFillings:
                #print('Related ', re.RelatedBuildingElement)
                if(re.RelatedBuildingElement.is_a()=='IfcWindow'):
                    windows.append(re.RelatedBuildingElement.id())
    return windows

def link_wall_window(ifcwalls):
    #link window and walls in plain python
    wallwindow=defaultdict(list)
    
    for wall in ifcwalls:
        wallwindow[wall.id()] = window_in_wall(wall)
        
    return wallwindow


class shadow_on_faces:
    """ simple container to hold computation results """
    def __init__(self,lfaces,lsun_dir):
        self._lfaces=lfaces
        self._lsun_dir=lsun_dir
        self._shadow_faces=[[] for i in range(len(self._lfaces))]
        self._durations_byfaces=[[]]
        

    def compute_shadow(self,exposed_building,min_area):
        for i,gf in enumerate(self._lfaces):
            # re computation of the face normal
            # shoudl be pointing outward
            srf = BRep_Tool().Surface(gf)
            plane = Geom_Plane.DownCast(srf)
            face_norm = plane.Axis().Direction()
            if(gf.Orientation()==1):
                face_norm.Reverse()
                
            for j,sun_dir in enumerate(self._lsun_dir):
                start=timer()
                shadow_face=shadow_caster_ext(sun_dir,exposed_building,gf,face_norm,1.e-3)
                print('     sun dir ',j,'/',len(self._lsun_dir))
                end=timer()
                self._shadow_faces[i].append(shadow_face)
                self._durations_byfaces[i].append(end-start)
                
        #print(' faces ',self._shadow_faces)
    
    def compute_area_and_ratio(self):
        gpp=GProp_GProps() 
        self._glass_area=0.0
        for gf in self._lfaces :
            brepgprop_SurfaceProperties(gf,gpp)
            self._glass_area+=gpp.Mass()
        
        self._shadow_area_vector=[]
        self._totalduration=[]
        for vector_idx in range(len(self._lsun_dir)):
            area_sum=0.0           
            for face_idx in range(len(self._lfaces)):
                brepgprop_SurfaceProperties(self._shadow_faces[face_idx][vector_idx],gpp)
                area_sum+=gpp.Mass()
                
            self._shadow_area_vector.append(area_sum)
            #self._totalduration.append( self._durations[face_idx])
            
        self._ratio_vector=[ a/self._glass_area for a in self._shadow_area_vector]
        #print(' shadow area vector ',self._shadow_area_vector)
        print(' ratio vector ',self._ratio_vector)
        
    def compute_area_and_ratio_byunion(self):
        """ 
        could be simpler in terms of code but rely on robustness of OCC to compute on more
        complex configurations 
        
        """
        totalface=fuse_listOfShape(self._lfaces)
        gpp=GProp_GProps() 
        brepgprop_SurfaceProperties(totalface,gpp)
        totalarea=gpp.Mass()
        
        ratio=[]
        for vector_idx in range(len(self._lsun_dir)):
            lfaces=[]
            for face_idx in range(len(self._lfaces)):
                f=self._shadow_faces[face_idx][vector_idx]
                if not f.IsNull():
                    lfaces.append(f)
            totalshadow=fuse_listOfShape(lfaces)
            brepgprop_SurfaceProperties(totalshadow,gpp)
            shadowarea=gpp.Mass()
            ratio.append(shadowarea/totalarea)
        
        #print(' shadow area vector ',self._shadow_area_vector)
        print(' ratio vector by union',ratio)
        
        
    def compute_complementary_face(self):
        cutter=BOPAlgo_BOP()
        
        self._complementary_faces=[[] for i in range(len(self._lfaces))]
        
        gpp=GProp_GProps() 
        
        #larea=[]
        
        for vector_idx in range(len(self._lsun_dir)):
            #area=0.0           
            for face_idx in range(len(self._lfaces)): 
                shadow_face=self._shadow_faces[face_idx][vector_idx]
                glass_face=self._lfaces[face_idx]
                
                if not shadow_face.IsNull():
                    cutter.Clear()
                    cutter.SetOperation(BOPAlgo_Operation.BOPAlgo_CUT)
                    cutter.AddArgument(glass_face)
                    cutter.AddTool(shadow_face)
                    cutter.SetFuzzyValue(1e-6)
                    cutter.Perform()
                    complementary=cutter.Shape()
                    #print(' cutter ',complementary)
                            
                else :
                    
                    complementary=glass_face
                    
                
                self._complementary_faces[face_idx].append(complementary)
                
               
class shadow_on_faces_byray:
    """ simple container to hold computation results """
    def __init__(self,lfaces,lsun_dir):
        self._lfaces=lfaces
        self._lsun_dir=lsun_dir
        self._shadow_tab=[[] for i in range(len(self._lfaces))]
        self._durations_byfaces=[[]]
        

    def compute_shadow(self,exposed_building,min_area,N):
        
        self._N=N
        for i,gf in enumerate(self._lfaces):
            # re computation of the face normal
            # shoudl be pointing outward
            srf = BRep_Tool().Surface(gf)
            plane = Geom_Plane.DownCast(srf)
            face_norm = plane.Axis().Direction()
            if(gf.Orientation()==1):
                face_norm.Reverse()
                
            for j,sun_dir in enumerate(self._lsun_dir):
                start=timer()
                #tab,lshape=shadow_caster_ray(sun_dir,exposed_building,gf,face_norm,N)
                tab=shadow_caster_ray(sun_dir,exposed_building,gf,face_norm,N)
                end=timer()
                self._shadow_tab[i].append(tab)
                #print(start,' ',end)
                self._durations_byfaces[i].append(end-start)
                
                
        #print(' faces ',self._shadow_faces)
    
    def compute_area_and_ratio(self):
        
        self._shadow_area_vector=[]
        self._totalduration=[]
        for vector_idx in range(len(self._lsun_dir)):
            area_sum=0.0           
            for face_idx in range(len(self._lfaces)):
                area_sum+=self._shadow_tab[face_idx][vector_idx].sum()
                
            self._shadow_area_vector.append(area_sum)
            #print(self._durations)
        
        
        #self._totalduration.append( self._durations[face_idx])
            
        self._ratio_vector=[ a/(self._N*self._N) for a in self._shadow_area_vector]
        #print(' shadow area vector ',self._shadow_area_vector)
        print(' ratio vector ray',self._ratio_vector)                
        
        

if __name__ == "__main__":

    # due to some bugs in ipython parsing
    __import__("logging").getLogger("parso.python.diff").setLevel("INFO")
    __import__("logging").getLogger("parso.cache").setLevel("INFO")
    __import__("logging").getLogger("asyncio").setLevel("INFO")

    window_tagname={
    'Ref':257076,
    'A1':257738,
    'A2':256772,
    'A3':257901,
    'B':256912,
    'C60':257017,
    'C45':266662
    }

    window_id_name={
    793:'Ref',
    840:'A1' ,
    566:'A2' ,
    857:'A3' ,
    757:'B'  ,
    776:'C60',
    998:'C45'
    }

    # Initialize a graphical display window (from ifcos)

    setting=ifcopenshell.geom.settings()
    setting.set(setting.USE_PYTHON_OPENCASCADE, True)

    ifc_file= ifcopenshell.open('data/Model_article_window.ifc')

    ifcwalls=ifc_file.by_type('IfcWall')
    ifcspaces=ifc_file.by_type('IfcSpace')
    ifcwindows=ifc_file.by_type('IfcWindow')
    ifcslabs=ifc_file.by_type('IfcSlab')
    ifcproxys=ifc_file.by_type('IfcBuildingElementProxy')
    """
    doors=ifc_file.by_type('IfcDoor')
    opening=ifc_file.by_type('IfcOpeningElement')
    storeys=ifc_file.by_type('IfcBuildingStorey')
    roof=ifc_file.by_type('IfcRoof')
    """
    tags=[w.Tag for w in ifcwindows]
    tags_ind=[ tags.index(str(t)) for t in window_tagname.values()]
    ifcwindows=[ifcwindows[indx] for indx in tags_ind]

    # partial building to compute external shell and exterior wall
    wall_shapes  = [create_shape(setting, x).geometry for x in ifcwalls if x.Representation is not None]
    space_shapes = [create_shape(setting, x).geometry for x in ifcspaces if x.Representation is not None]
    core_shapes  = wall_shapes+space_shapes
    core_solids  = shapes_as_solids(core_shapes)
        
    # complete building to compute shadow on
    ifcextension= []+ifcslabs+ifcproxys
    extension_shapes = [create_shape(setting, x).geometry for x in ifcextension if x.Representation is not None]
    extension_solids =  shapes_as_solids(extension_shapes)

    
    building_shapes= core_shapes + extension_shapes
    building_solids= core_solids + extension_solids
    exposed_building = fuse_listOfShape(building_solids)

    external_shell= get_external_shell(core_solids)
    
    windows_by_wall_id = dict()
    for wall in ifcwalls:
        windows_by_wall_id[wall.id()] = window_in_wall(wall)

    # will only contain wall id that considered as exterior wall of the building
    # if id is not in the keys, the wall could be considered as interior
    normal_by_exterior_wall_id = dict()
    for (w_id,ws) in zip(windows_by_wall_id.keys(),wall_shapes):
        normal_by_exterior_wall_id[w_id]=exterior_wall_normal(ws,external_shell)
    
    
    glassface_bywindowid=defaultdict(list)
    for w_id in windows_by_wall_id.keys():
        if (w_id in normal_by_exterior_wall_id.keys()):
            wall_norm=normal_by_exterior_wall_id[w_id]
        else:
            # window in interior wall
            continue
            
        for win_id in windows_by_wall_id[w_id]:
        
            windowshape=create_shape(setting, ifc_file.by_guid(win_id)).geometry
            gfaces=biggestface_along_vector(windowshape,wall_norm)
            glassface_bywindowid[win_id].extend(gfaces)
    
    
    npos=5
    h_angles=np.arange(0,360.,360./npos)
    v_angles=[60.]#,65.,70.,75.,80.,85.]


    x=-np.cos(np.deg2rad(h_angles))
    y=-np.sin(np.deg2rad(h_angles))
    z=-np.sin(np.deg2rad(v_angles))

    lparams=[(v,h) for (v,h) in itertools.product(v_angles,h_angles)]
    vvalues=[ v[0] for v in lparams]
    hvalues=[ v[1] for v in lparams]

    l_sun_dir=[gp_Dir(xi,yi,zi) for (zi,(xi,yi)) in itertools.product(z,zip(x,y))]

    N=[4,6,8,10,15,20,30,40,50]               
            
    lsof=[]

    for (k,win_id) in enumerate(glassface_bywindowid.keys()):
        lglassfaces=glassface_bywindowid[win_id]
        print(' window id ', win_id)
        sof=shadow_on_faces(lglassfaces,l_sun_dir)
        sof.compute_shadow(exposed_building,1e-3)
        sof.compute_area_and_ratio()
        sof.compute_complementary_face()
        lsof.append(sof)




    """
    llsofr=[[] for nray in N]    
    for i,nray in enumerate(N):
        
        for (k,win_id) in enumerate(glassface_bywindowid.keys()):
            lglassfaces=glassface_bywindowid[win_id]
            sofr = shadow_on_faces_byray(lglassfaces,l_sun_dir)
            sofr.compute_shadow(exposed_building,1e-3,nray)
            sofr.compute_area_and_ratio()
            llsofr[i].append((nray,sofr))
    """
    # Build a dataframe with all the results
    frames=[]
    for sof,id in zip(lsof,glassface_bywindowid.keys()):
        name=window_id_name[id]
        durations=np.array(sof._durations_byfaces).mean(axis=0)
        frames.append(pd.DataFrame(zip(vvalues,hvalues,
                                    itertools.repeat(name),
                                    itertools.repeat(0),
                                    sof._ratio_vector,
                                    durations)))
    """
    for i,(lsofr,n) in enumerate(zip(llsofr,N)):
        
        for (_,sof),id in zip(lsofr,glassface_bywindowid.keys()):
            name=window_id_name[id]
            durations=np.array(sof._durations_byfaces).mean(axis=0)
            frames.append(pd.DataFrame(zip(vvalues,hvalues,
                                    itertools.repeat(name),
                                    itertools.repeat(n),
                                    sof._ratio_vector,
                                    durations)))
    """  
    res=pd.concat(frames) 
    res.columns=['v_angle','h_angle','name','Nray','shad_ratio','duration']
    res.to_csv('Results_extrusiononly.csv')


    """
    # to export some shadow in stl files
    # work well for unique shadow by window (no staking)
    lshape=[s for sof in lsof for s in sof._shadow_faces]
    lshape=list(itertools.chain(*lshape))
    write_stl_file(shadowshape, 'shadow.stl')

    lshape=[s for sof in lsof for s in sof._complementary_faces]
    lshape=list(itertools.chain(*lshape))
    shadowshape=fuse_listOfShape(lshape)
    write_stl_file(shadowshape, 'shadow_compl.stl')
    """


    """

    def rgb_color(r, g, b):
        return Quantity_Color(r, g, b, Quantity_TOC_RGB)
    x=50/256
    gray=rgb_color(x, x, x)

    display, start_display, add_menu, add_function_to_menu = init_display()

    [display.DisplayShape(s,color=gray,transparency=0.9) for s in building_shapes]
    for sof in lsof:
        [display.DisplayShape(s,transparency=0.1,color='BLACK') for s in sof._shadow_faces]
        [display.DisplayShape(s,transparency=0.1,color='YELLOW') for s in sof._complementary_faces]
        
    [display.DisplayShape(s,color='RED',transparency=0.0) for s in lext2] 
    #[display.DisplayShape(s,color='GREEN',transparency=0.1) for s in lshells] 

    [display.DisplayShape(s,color='BLUE',transparency=0.1) for s in lf] 
    #[display.DisplayShape(s,color='RED',transparency=0.0) for s in lcyl]

    #s=lsof[0]._complementary_faces[0]

    display.FitAll()
    #ifcopenshell.geom.utils.main_loop()
    start_display()
    

    """

