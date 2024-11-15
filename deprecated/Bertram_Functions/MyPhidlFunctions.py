
import numpy as np
from phidl import LayerSet
from phidl import quickplot as qp
from phidl import Path, CrossSection, Device
import phidl.path as pp
import phidl.geometry as pg
import phidl


def WaveGuideMaker(SegWidth, SegSpacing, SegLength, SegPorts):
    #First create the path
    P = Path()
    P.append(pp.straight(length = SegLength))

    #Then create the cross section
    CS = CrossSection()
    CS.add(width=SegWidth, offset = 0, layer = 1, name = 'CenterConductor')
    CS.add(width=SegWidth + 2*SegSpacing, offset = 0, layer= 1, name = 'Etch', ports= [SegPorts[0],SegPorts[1]])

    #Finally we extrude
    Device = P.extrude(CS)
    
    return Device, CS


def SegmentMaker(SegWidth, SegLength, SegLayer, SegName, SegPorts):
    #First create the path
    P = Path()
    P.append(pp.straight(length = SegLength))

    #Then create the cross section
    CS = CrossSection()
    CS.add(width=SegWidth, offset = 0, layer= SegLayer, name = SegName, ports= [SegPorts[0],SegPorts[1]])

    #Finally we extrude
    Device = P.extrude(CS)
    
    return Device, CS

def TransitionSegmentMaker(CrossSection1,CrossSection2,SegLength):
    CStrans = pp.transition(cross_section1 = CrossSection1,
                            cross_section2 = CrossSection2,
                            width_type = 'linear')
    P = pp.straight(length = SegLength)
    Device = P.extrude(CStrans)

    return Device

def FeedLineMakerOLD(FeedLineLength,FeedLineWidth,BondPadLength,BondPadWidth,SpacingFraction,BondPadToFeedLineTransLength):
    #FeedLine and etched area on top and below feedline. Needs a heal and a subtract in beamer.
    FeedLineSeg, FeedLineCS = WaveGuideMaker(FeedLineWidth, SpacingFraction*FeedLineWidth, FeedLineLength, SegPorts = ['FL_in','FL_out'])
    BondPadLeftSeg, BondPadLeftCS = WaveGuideMaker(BondPadWidth, SpacingFraction*BondPadWidth, BondPadLength, SegPorts = ['BPL_in','BPL_out'])
    BondPadRightSeg, BondPadRightCS = WaveGuideMaker(BondPadWidth, SpacingFraction*BondPadWidth, BondPadLength, SegPorts = ['BPR_in','BPR_out'])

    #Transition from bondpad to feedline
    BondPadLeftToFeedLineTransSeg = TransitionSegmentMaker(BondPadLeftCS,FeedLineCS,BondPadToFeedLineTransLength)
    FeedLineToBondPadTightTransSeg = TransitionSegmentMaker(FeedLineCS,BondPadRightCS,BondPadToFeedLineTransLength)

    #Etched area at the end of the feedlines
    EtchAreaLeftSeg,EtchAreaLeftCS = SegmentMaker(BondPadWidth + 2*BondPadWidth*SpacingFraction, BondPadWidth*SpacingFraction, 1, 'Etch', [None,'EAL_out'])
    EtchAreaRightSeg,EtchAreaRightCS = SegmentMaker(BondPadWidth + 2*BondPadWidth*SpacingFraction, BondPadWidth*SpacingFraction, 1, 'Etch', ['EAR_in',None])

    D = Device()

    FeedLine = D << FeedLineSeg
    BondPadLeft = D << BondPadLeftSeg
    BondPadRight = D << BondPadRightSeg
    BondPadLeftToFeedLineTrans = D << BondPadLeftToFeedLineTransSeg
    FeedLineToBondPadTightTrans = D << FeedLineToBondPadTightTransSeg
    EtchAreaLeft = D << EtchAreaLeftSeg
    EtchAreaRight = D << EtchAreaRightSeg


    BondPadLeft.connect('BPL_in',EtchAreaLeft.ports['EAL_out'])
    BondPadLeftToFeedLineTrans.connect('FL_in',BondPadLeft.ports['BPL_out'])
    FeedLine.connect('FL_in', BondPadLeftToFeedLineTrans.ports['BPL_out'])
    FeedLineToBondPadTightTrans.connect('BPR_in', FeedLine.ports['FL_out'])
    BondPadRight.connect('BPR_in', FeedLineToBondPadTightTrans.ports['FL_out'])
    EtchAreaRight.connect('EAR_in', BondPadRight.ports['BPR_out'])

    return BondPadRight,EtchAreaRight


def FeedLineMaker(FeedLineWidth = None,
                  FeedLineSpacing = None,
                  FeedLineLength = None,
                  TransitionLength = None,
                  BondPadLength = None,
                  BondPadWidth = None):

    #### Sample Feed Line Values ####
    if FeedLineWidth == None:
        FeedLineWidth = 95
    if FeedLineSpacing == None:
        FeedLineSpacing = 16
    if FeedLineLength == None:
        FeedLineLength = 6000
    if TransitionLength == None:
        TransitionLength = 200
    if BondPadLength == None:
        BondPadLength = 300
    if BondPadWidth == None:
        BondPadWidth = 300
   
    ####Metal Area #####

    FeedLine_device = Device()
    FeedLine_rect = pg.taper(length = FeedLineLength, width1 = FeedLineWidth, width2 = FeedLineWidth, port = None, layer = 0)

    TransitionRight_rect = pg.taper(length = TransitionLength, width1 = BondPadWidth, width2 = FeedLineWidth, port = None, layer = 0)
    TransitionLeft_rect = pg.taper(length = TransitionLength, width1 = FeedLineWidth, width2 = BondPadWidth, port = None, layer = 0)

    BondPadRight_rect = pg.taper(length = BondPadLength, width1 = BondPadWidth, width2 = BondPadWidth, port = None, layer = 0)
    BondPadLeft_rect = pg.taper(length = BondPadLength, width1 = BondPadWidth, width2 = BondPadWidth, port = None, layer = 0)

    FeedLine_rect_ref = FeedLine_device << FeedLine_rect
    TransitionRight_rect_ref = FeedLine_device << TransitionRight_rect
    TransitionLeft_rect_ref = FeedLine_device << TransitionLeft_rect
    BondPadRight_rect_ref = FeedLine_device << BondPadRight_rect
    BondPadLeft_rect_ref = FeedLine_device << BondPadLeft_rect

    TransitionLeft_rect_ref.connect(2,BondPadLeft_rect_ref.ports[1])
    FeedLine_rect_ref.connect(2,TransitionLeft_rect_ref.ports[1])
    TransitionRight_rect_ref.connect(2,FeedLine_rect_ref.ports[1])
    BondPadRight_rect_ref.connect(2,TransitionRight_rect_ref.ports[1])

    FeedLine_device_CenterCoord = FeedLine_device.center
    FeedLine_device.move(0-FeedLine_device_CenterCoord)

    ####Etched Area#####

    FeedLineEtch_device = Device()

    SpacingFraction = FeedLineSpacing/FeedLineWidth
    BondPadSpacing = BondPadWidth*SpacingFraction

    FeedLineEtch_rect = pg.taper(length = FeedLineLength, width1 = FeedLineWidth+FeedLineSpacing*2, width2 = FeedLineWidth+FeedLineSpacing*2, port = None, layer = 0)

    TransitionRightEtch_rect = pg.taper(length = TransitionLength, width1 = BondPadWidth+BondPadSpacing*2 , width2 = FeedLineWidth+FeedLineSpacing*2, port = None, layer = 0)
    TransitionLeftEtch_rect = pg.taper(length = TransitionLength, width1 = FeedLineWidth+FeedLineSpacing*2, width2 = BondPadWidth+BondPadSpacing*2, port = None, layer = 0)

    BondPadRightEtch_rect = pg.taper(length = BondPadLength, width1 = BondPadWidth+BondPadSpacing*2, width2 = BondPadWidth+BondPadSpacing*2, port = None, layer = 0)
    BondPadLeftEtch_rect = pg.taper(length = BondPadLength, width1 = BondPadWidth+BondPadSpacing*2, width2 = BondPadWidth+BondPadSpacing*2, port = None, layer = 0)

    EndEtchRight_rect = pg.taper(length = BondPadSpacing, width1 = BondPadWidth+BondPadSpacing*2, width2 = BondPadWidth+BondPadSpacing*2, port = None, layer = 0)
    EndEtchLeft_rect = pg.taper(length = BondPadSpacing, width1 = BondPadWidth+BondPadSpacing*2, width2 = BondPadWidth+BondPadSpacing*2, port = None, layer = 0)


    FeedLineEtch_rect_ref = FeedLineEtch_device << FeedLineEtch_rect
    TransitionRightEtch_rect_ref = FeedLineEtch_device << TransitionRightEtch_rect
    TransitionLeftEtch_rect_ref = FeedLineEtch_device << TransitionLeftEtch_rect
    BondPadRightEtch_rect_ref = FeedLineEtch_device << BondPadRightEtch_rect
    BondPadLeftEtch_rect_ref = FeedLineEtch_device << BondPadLeftEtch_rect
    EndEtchRight_rect_ref = FeedLineEtch_device << EndEtchRight_rect
    EndEtchLeft_rect_ref = FeedLineEtch_device << EndEtchLeft_rect

    BondPadLeftEtch_rect_ref.connect(2,EndEtchLeft_rect_ref.ports[1])
    TransitionLeftEtch_rect_ref.connect(2,BondPadLeftEtch_rect_ref.ports[1])
    FeedLineEtch_rect_ref.connect(2,TransitionLeftEtch_rect_ref.ports[1])
    TransitionRightEtch_rect_ref.connect(2,FeedLineEtch_rect_ref.ports[1])
    BondPadRightEtch_rect_ref.connect(2,TransitionRightEtch_rect_ref.ports[1])
    EndEtchRight_rect_ref.connect(2,BondPadRightEtch_rect_ref.ports[1])

    FeedLineEtch_device_CenterCoord = FeedLineEtch_device.center
    FeedLineEtch_device.move(0-FeedLineEtch_device_CenterCoord)

    #### Subtract ####

    EtchedArea = pg.xor_diff(A = FeedLineEtch_device,B = FeedLine_device, precision=1e-6)

    return EtchedArea

def ArrayRotator(array):
    index = (np.sum(array,axis=1).tolist()).index(max(np.sum(array,axis=1)))
    return np.roll(array,-index,axis = 0)


def FourProbeMaker(BarLength, BarWidth, BondPadLength, BondPadWidth, ProbeLength, ProbeWidth, Spacing):
    
    FourProbeDevice = Device()

    BarDevice = Device()
    Bar_temp = pg.rectangle([BarLength,BarWidth], layer = 0)
    Bar = BarDevice << Bar_temp

    BarDevice.add_port(name = 'Bar_Left', midpoint = [0,BarWidth/2], width = BarWidth, orientation = 180)
    BarDevice.add_port(name = 'Bar_Right', midpoint = [BarLength,BarWidth/2], width = BarWidth, orientation = 0)

    BarDevice.add_port(name = 'Probe_Left', midpoint = [BarLength/3,BarWidth], width = ProbeWidth, orientation = 90)
    BarDevice.add_port(name = 'Probe_Right', midpoint = [(BarLength*2)/3,BarWidth], width = ProbeWidth, orientation = 90)

    ProbeDevice = Device()
    Probe_temp = pg.rectangle([ProbeLength,ProbeWidth], layer = 0)
    Probe = ProbeDevice << Probe_temp

    ProbeDevice.add_port(name = 'Probe_Bottom', midpoint = [0,ProbeWidth/2], width = ProbeWidth, orientation = 180)
    ProbeDevice.add_port(name = 'Probe_Top', midpoint = [ProbeLength,ProbeWidth/2], width = ProbeWidth, orientation = 0)

    BondPadDevice = Device()
    BondPad_temp = pg.rectangle([BondPadLength,BondPadWidth], layer = 0)
    BondPad = BondPadDevice << BondPad_temp

    BondPadDevice.add_port(name = 'BondPad_Port', midpoint=[0,BondPadWidth/2], width = BondPadWidth, orientation = 180)

    ADevice1 = Device()
    BarRef = ADevice1 << BarDevice
    ProbeLeftRef = ADevice1 << ProbeDevice
    ProbeRightRef = ADevice1 << ProbeDevice
    BondPadLeftRef = ADevice1 << BondPadDevice
    BondPadProbeLeftRef = ADevice1 << BondPadDevice
    BondPadProbeRighttRef = ADevice1 << BondPadDevice
    BondPadRightRef = ADevice1 << BondPadDevice

    BondPadLeftRef.connect(port = 'BondPad_Port', destination = BarRef.ports['Bar_Left'])
    BondPadRightRef.connect(port = 'BondPad_Port', destination = BarRef.ports['Bar_Right'])
    ProbeLeftRef.connect(port = 'Probe_Bottom', destination = BarRef.ports['Probe_Left'])
    ProbeRightRef.connect(port = 'Probe_Bottom', destination = BarRef.ports['Probe_Right'])
    BondPadProbeLeftRef.connect(port = 'BondPad_Port', destination = ProbeLeftRef.ports['Probe_Top'])
    BondPadProbeRighttRef.connect(port = 'BondPad_Port', destination = ProbeRightRef.ports['Probe_Top'])

    FourProbeArea = pg.union(ADevice1,layer = 0)
    FourProbeDevice << FourProbeArea

    SpacingArray = np.array([[1,1],
                            [1,-1],
                            [-1,-1],
                            [-1,1]]) * Spacing

    PolygonArray = [ArrayRotator(polygon) + SpacingArray  for polygon in ADevice1.get_polygons()]

    ADevice2 = Device()
    ADevice2.add_polygon(PolygonArray,layer=1)
    EtchArea = pg.union(ADevice2, layer = 1)

    FourProbeEtch = FourProbeDevice << EtchArea

    return FourProbeDevice



def SchusterResonator(NumberOfBends = None
                      ,HorizontalBarLength = None
                      ,VerticalBarLength = None
                      ,EndTailLength = None
                      ,ResonatorWidth = None
                      ,CouplingCapacitorWidth = None
                      ,CouplingCapacitorTopBarLength = None
                      ,CouplingCapacitorSideWallsLength = None
                      ,Spacing = None
                      ,NanowireArea = None):

    #Resonator
    if NumberOfBends is None:
        NumberOfBends = 16
    if HorizontalBarLength is None:
        HorizontalBarLength = 40
    if VerticalBarLength is None:
        VerticalBarLength = 5
    if EndTailLength is None:
        EndTailLength = 30
    if ResonatorWidth is None:
        ResonatorWidth = 1

    #Coupling Capacitor
    if CouplingCapacitorWidth is None:
        CouplingCapacitorWidth = 4
    if CouplingCapacitorTopBarLength is None:
        CouplingCapacitorTopBarLength = 60
    if CouplingCapacitorSideWallsLength is None:
        CouplingCapacitorSideWallsLength = 60

    #Spacing for etch area
    if Spacing is None:
        Spacing = 30

    #AreaForNanowire
    if NanowireArea is None:
        NanowireArea = 0

    #Making the resonator
    P_Resonator = Path()

    P_Resonator.append(pp.straight(length = VerticalBarLength, num_pts = 2))
    P_Resonator.end_angle += 90
    P_Resonator.append(pp.straight(length = HorizontalBarLength/2, num_pts = 2))
    P_Resonator.end_angle += -90
    P_Resonator.append(pp.straight(length=VerticalBarLength, num_pts = 2))
    P_Resonator.end_angle += -90

    for i in range(NumberOfBends):
        P_Resonator.append(pp.straight(length=HorizontalBarLength, num_pts = 2))
        P_Resonator.end_angle += 90*((-1)**(i%2))
        P_Resonator.append(pp.straight(length=VerticalBarLength, num_pts = 2))
        P_Resonator.end_angle += 90*((-1)**(i%2))

    P_Resonator.append(pp.straight(length=HorizontalBarLength/2, num_pts = 2))
    P_Resonator.end_angle += 90-180*(NumberOfBends%2)
    P_Resonator.append(pp.straight(length=EndTailLength, num_pts = 2))

    P_Resonator.rotate(-90)

    X_Resonator = CrossSection()

    X_Resonator.add(width = ResonatorWidth, offset = 0, layer = 0)
    Resonator_Device = P_Resonator.extrude(X_Resonator)

    Resonator_Device.add_port(name='Top', midpoint=(0, 0), width=ResonatorWidth, orientation=90, port=None)


    #Making the Coupling capacitor

    P_CouplingCapacitor = Path()

    P_CouplingCapacitor.append(pp.straight(length = CouplingCapacitorSideWallsLength, num_pts = 2))
    P_CouplingCapacitor.end_angle += 90
    P_CouplingCapacitor.append(pp.straight(length = CouplingCapacitorTopBarLength, num_pts = 2))
    P_CouplingCapacitor.end_angle += 90
    P_CouplingCapacitor.append(pp.straight(length = CouplingCapacitorSideWallsLength, num_pts = 2))

    P_CouplingCapacitor.rotate(90)

    X_CouplingCapacitor = CrossSection()
    X_CouplingCapacitor.add(width = CouplingCapacitorWidth, offset = 0, layer = 0)
    CouplingCapacitor_device = P_CouplingCapacitor.extrude(X_CouplingCapacitor)

    CouplingCapacitor_device.add_port(name='Bottom', midpoint=(-CouplingCapacitorTopBarLength/2, CouplingCapacitorSideWallsLength-CouplingCapacitorWidth/2), width=ResonatorWidth, orientation=270, port=None)


    SchusterResonatorDevice = Device()
    CouplingCapacitor_Ref = SchusterResonatorDevice << CouplingCapacitor_device
    Resonator_Device_Ref = SchusterResonatorDevice << Resonator_Device

    Resonator_Device_Ref.connect(port = 'Top', destination = CouplingCapacitor_Ref.ports['Bottom'])

    #EtchArea
    EtchAreaBoundary_Device = Device()

    BoundingBox = SchusterResonatorDevice.bbox
    
    BoundingBox[0,0] += -Spacing
    BoundingBox[1,0] += Spacing
    BoundingBox[1,1] += Spacing
    print()

    BoundingBox_Xpts = (BoundingBox[0,0],BoundingBox[0,0],BoundingBox[1,0],BoundingBox[1,0])
    BoundingBox_Ypts = (BoundingBox[0,1] - NanowireArea,BoundingBox[1,1],BoundingBox[1,1],BoundingBox[0,1] - NanowireArea)

    EtchAreaBoundary_Device.add_polygon([BoundingBox_Xpts,BoundingBox_Ypts],layer=0)

    EtchedArea = pg.xor_diff(A = EtchAreaBoundary_Device,B = SchusterResonatorDevice, precision=1e-6)

    return EtchedArea, SchusterResonatorDevice, EtchAreaBoundary_Device, P_Resonator.length()


def ResonatorMover(Chip_ref,Resonator_ref,MoveVectorFromCenter):

    ResonatorDevice1_BoundingBox = Resonator_ref.bbox
    ResonatorDevice1_TopCoord = [np.average(ResonatorDevice1_BoundingBox,axis=0)[0],ResonatorDevice1_BoundingBox[1][1]]
    MoveVector = -np.array(ResonatorDevice1_TopCoord)+np.array(Chip_ref.center)+np.array(MoveVectorFromCenter)

    return MoveVector

def MarkerFieldMaker(MarkerLength,MarkerWidth,MarkerFieldWidth,MarkerFieldLength):
    MarkerDevice = Device()

    MarkerTR = pg.cross(length = MarkerLength, width = MarkerWidth, layer = 0)
    MarkerBR = pg.cross(length = MarkerLength, width = MarkerWidth, layer = 0)
    MarkerBL = pg.cross(length = MarkerLength, width = MarkerWidth, layer = 0)
    MarkerTL = pg.cross(length = MarkerLength, width = MarkerWidth, layer = 0)

    MarkerTR.move([MarkerFieldWidth,MarkerFieldLength])
    MarkerBR.move([MarkerFieldWidth,0])
    MarkerBL.move([0,0])
    MarkerTL.move([0,MarkerFieldLength])


    MarkerTR_ref = MarkerDevice << MarkerTR
    MarkerBR_ref = MarkerDevice << MarkerBR
    MarkerBL_ref = MarkerDevice << MarkerBL
    MarkerTL_ref = MarkerDevice << MarkerTL

    return MarkerDevice

def MarkerArrayIdentifier_func(ArraySize,textScale,spacing,crossLength,crossWidth):

    D = Device()
    for i in range(ArraySize[0]):
        for j in range(ArraySize[1]):
            Cross_Device = pg.cross(length=crossLength,width=crossWidth, layer = 0)
            Cross_Device.move(np.multiply(np.array(spacing),np.array([i,j])))
            D << Cross_Device

    for i in range(ArraySize[0]):
        for j in range(ArraySize[1]):
            Indentifier_Device = pg.text(str([i,j]),size=textScale,layer = 0)
            Indentifier_Device.move((np.multiply(np.array(spacing),np.array([i,j])))+np.array([crossWidth,crossWidth])*2)
            D << Indentifier_Device

    return D

def MarkerArray_func(ArraySize,spacing,crossLength,crossWidth):

    D = Device()
    for i in range(ArraySize[0]):
        for j in range(ArraySize[1]):
            Cross_Device = pg.cross(length=crossLength,width=crossWidth, layer = 0)
            Cross_Device.move(np.multiply(np.array(spacing),np.array([i,j])))
            D << Cross_Device

    return D

def InnerOuterMarkerField(spacingInner,crossLengthInner,crossWidthInner,spacingOuter,crossLengthOuter,crossWidthOuter):
    D = Device()
    for i in range(2):
        for j in range(2):
            CrossInner_poly = pg.cross(length=crossLengthInner,width=crossWidthInner, layer = 0)
            if [i,j] == [0,0]:
                CrossInner_poly.move(np.array([(spacingOuter/2)-(spacingInner/2),(spacingOuter/2)-(spacingInner/2)]))
            elif [i,j] == [0,1]:
                CrossInner_poly.move(np.array([(spacingOuter/2)-(spacingInner/2),(spacingOuter/2)+(spacingInner/2)]))
            elif [i,j] == [1,1]:
                CrossInner_poly.move(np.array([(spacingOuter/2)+(spacingInner/2),(spacingOuter/2)+(spacingInner/2)]))
            else:
                CrossInner_poly.move(np.array([(spacingOuter/2)+(spacingInner/2),(spacingOuter/2)-(spacingInner/2)]))

            D << CrossInner_poly

    for i in range(2):
        for j in range(2):
            Cross_Device = pg.cross(length=crossLengthOuter,width=crossWidthOuter, layer = 0)
            Cross_Device.move(np.multiply(np.array(spacingOuter),np.array([i,j])))
            D << Cross_Device

    return D

def FeedLineResonator(ResonatorLength = None,
                      ResonatorWidth  = None,
                      FeedLineWidth   = None,
                      FeedLineSpacing = None,
                      TotalFeedLineLength  = None,
                      BondPadLength = None,
                      BondPadWidth = None,
                      TransitionLength = None):

    if ResonatorLength == None:
        ResonatorLength = 1500
    if ResonatorWidth == None:
        ResonatorWidth  = 2
    if FeedLineWidth == None:
        FeedLineWidth   = 85
    if FeedLineSpacing == None:
        FeedLineSpacing = 2
    if TotalFeedLineLength == None:
        TotalFeedLineLength  = 7000
    if BondPadLength == None:
        BondPadLength = 300
    if BondPadWidth == None:
        BondPadWidth = 300
    if TransitionLength == None:
        TransitionLength = 200

    #Inner pin of the resonator
    BondPadRight_rect = pg.taper(length = BondPadLength, width1 = BondPadWidth, width2 = BondPadWidth, port = None, layer = 0)
    BondPadLeft_rect = pg.taper(length = BondPadLength, width1 = BondPadWidth, width2 = BondPadWidth, port = None, layer = 0)

    TransitionLeft_rect = pg.taper(length = TransitionLength, width1 = BondPadWidth, width2 = FeedLineWidth, port = None, layer = 0)
    TransitionRight_rect = pg.taper(length = TransitionLength, width1 = FeedLineWidth, width2 = BondPadWidth, port = None, layer = 0)

    FeedLine_rect = pg.taper(length = (TotalFeedLineLength - ResonatorLength - 2*TransitionLength - 2*BondPadLength)/2, width1 = FeedLineWidth, width2 = FeedLineWidth, port = None, layer = 0)
    Resonator_rect = pg.taper(length = ResonatorLength, width1 = ResonatorWidth, width2 = ResonatorWidth, port = None, layer = 0)


    FeedLineResonatorInnerPin_device = Device()

    BondPadRight_ref = FeedLineResonatorInnerPin_device << BondPadRight_rect
    BondPadLeft_ref = FeedLineResonatorInnerPin_device << BondPadLeft_rect
    TransitionRight_ref = FeedLineResonatorInnerPin_device << TransitionRight_rect
    TransitionLeft_ref = FeedLineResonatorInnerPin_device << TransitionLeft_rect
    FeedLineRight_ref = FeedLineResonatorInnerPin_device << FeedLine_rect
    FeedLineLeft_ref = FeedLineResonatorInnerPin_device << FeedLine_rect
    Resonator_ref = FeedLineResonatorInnerPin_device << Resonator_rect

    TransitionLeft_ref.connect(1,BondPadLeft_ref.ports[2])
    FeedLineLeft_ref.connect(1,TransitionLeft_ref.ports[2])
    Resonator_ref.connect(1,FeedLineLeft_ref.ports[2])
    FeedLineRight_ref.connect(1,Resonator_ref.ports[2])
    TransitionRight_ref.connect(1,FeedLineRight_ref.ports[2])
    BondPadRight_ref.connect(1,TransitionRight_ref.ports[2])

    #Spacing shape

    SpacingFraction = FeedLineSpacing/FeedLineWidth
    BondPadSpacing = SpacingFraction*BondPadWidth


    EndPaddingRight_rect = pg.taper(length = BondPadSpacing, width1 = BondPadWidth+2*BondPadSpacing, width2 = BondPadWidth+2*BondPadSpacing, port = None, layer = 0)
    EndPaddingLeft_rect = pg.taper(length = BondPadSpacing, width1 = BondPadWidth+2*BondPadSpacing, width2 = BondPadWidth+2*BondPadSpacing, port = None, layer = 0)
    BondPadSpacingRight_rect = pg.taper(length = BondPadLength, width1 = BondPadWidth+2*BondPadSpacing, width2 = BondPadWidth+2*BondPadSpacing, port = None, layer = 0)
    BondPadSpacingLeft_rect = pg.taper(length = BondPadLength, width1 = BondPadWidth+2*BondPadSpacing, width2 = BondPadWidth+2*BondPadSpacing, port = None, layer = 0)

    TransitionSpacingLeft_rect = pg.taper(length = TransitionLength, width1 = BondPadWidth+2*BondPadSpacing, width2 = FeedLineWidth+2*FeedLineSpacing, port = None, layer = 0)
    TransitionSpacingRight_rect = pg.taper(length = TransitionLength, width1 = FeedLineWidth+2*FeedLineSpacing, width2 = BondPadWidth+2*BondPadSpacing, port = None, layer = 0)

    FeedLineSpacing_rect = pg.taper(length = TotalFeedLineLength - 2*TransitionLength - 2*BondPadLength, width1 = FeedLineWidth+2*FeedLineSpacing, width2 = FeedLineWidth+2*FeedLineSpacing, port = None, layer = 0)

    FeedLineResonatorSpacing_device = Device()

    EndPaddingRight_ref = FeedLineResonatorSpacing_device << EndPaddingRight_rect
    EndPaddingLeft_ref = FeedLineResonatorSpacing_device << EndPaddingLeft_rect
    BondPadSpacingRight_ref = FeedLineResonatorSpacing_device << BondPadSpacingRight_rect
    BondPadSpacingLeft_ref = FeedLineResonatorSpacing_device << BondPadSpacingLeft_rect
    TransitionSpacingLeft_ref = FeedLineResonatorSpacing_device << TransitionSpacingLeft_rect
    TransitionSpacingRight_ref = FeedLineResonatorSpacing_device << TransitionSpacingRight_rect
    FeedLineSpacing_ref = FeedLineResonatorSpacing_device << FeedLineSpacing_rect


    BondPadSpacingLeft_ref.connect(1,EndPaddingLeft_ref.ports[2])
    TransitionSpacingLeft_ref.connect(1,BondPadSpacingLeft_ref.ports[2])
    FeedLineSpacing_ref.connect(1,TransitionSpacingLeft_ref.ports[2])
    TransitionSpacingRight_ref.connect(1,FeedLineSpacing_ref.ports[2])
    BondPadSpacingRight_ref.connect(1,TransitionSpacingRight_ref.ports[2])
    EndPaddingRight_ref.connect(1,BondPadSpacingRight_ref.ports[2])

    FeedLineResonatorInnerPin_device.movex(BondPadSpacing)
    EtchedArea = pg.xor_diff(A = FeedLineResonatorSpacing_device,B = FeedLineResonatorInnerPin_device, precision=1e-6)

    return EtchedArea, FeedLineResonatorSpacing_device, FeedLineResonatorInnerPin_device




def FeedlineResonatorForSim(ResonatorLength = None,
                      ResonatorWidth = None,
                      FeedLineWidth = None,
                      FeedLineLength = None,
                      FeedLineSpacing = None,
                      ChipHeight = None):
    
    if ResonatorLength == None:
        ResonatorLength = 1500
    if ResonatorWidth == None:
        ResonatorWidth  = 2
    if FeedLineWidth == None:
        FeedLineWidth   = 85
    if FeedLineLength == None:
        FeedLineLength  = 2000
    if FeedLineSpacing == None:
        FeedLineSpacing = 2
    if ChipHeight == None:
        ChipHeight      = 2000

    FeedLineResonatorInnerPin_device = Device()
    FeedLine_rect = pg.taper(length = FeedLineLength, width1 = FeedLineWidth, width2 = FeedLineWidth, port = None, layer = 0)
    Resonator_rect = pg.taper(length = ResonatorLength, width1 = ResonatorWidth, width2 = ResonatorWidth, port = None, layer = 0)

    FeedLineRight_ref = FeedLineResonatorInnerPin_device << FeedLine_rect
    FeedLineLeft_ref = FeedLineResonatorInnerPin_device << FeedLine_rect
    Resonator_ref = FeedLineResonatorInnerPin_device << Resonator_rect

    Resonator_ref.connect(2,FeedLineLeft_ref.ports[1])
    FeedLineRight_ref.connect(2,Resonator_ref.ports[1])

    FeedLineResonator_Device = Device()

    Spacing_rect = pg.taper(length = 2*FeedLineLength+ResonatorLength+2*FeedLineSpacing, width1 = FeedLineWidth+2*FeedLineSpacing, width2 = FeedLineWidth+2*FeedLineSpacing, port = None, layer = 0)

    Chip_rect = pg.rectangle(size=(2*FeedLineLength+ResonatorLength+2*FeedLineSpacing, ChipHeight))


    FeedLineResonatorInnerPin_ref = FeedLineResonator_Device << FeedLineResonatorInnerPin_device
    Spacing_ref = FeedLineResonator_Device << Spacing_rect
    Chip_ref = FeedLineResonator_Device << Chip_rect

    FeedLineResonatorInnerPin_ref.move(Chip_ref.center - FeedLineResonatorInnerPin_ref.center)
    Spacing_ref.move(Chip_ref.center - Spacing_ref.center)


    EtchedArea = pg.xor_diff(A = Spacing_ref,B = Chip_ref, precision=1e-6)

    FinalDevice = Device()
    EtchedArea = FinalDevice << EtchedArea
    FeedLineResonatorInnerPin_ref = FinalDevice << FeedLineResonatorInnerPin_device
    FeedLineResonatorInnerPin_ref.move(EtchedArea.center-FeedLineResonatorInnerPin_ref.center)

    return FinalDevice
