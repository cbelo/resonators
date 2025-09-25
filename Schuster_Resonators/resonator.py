import numpy as np
from phidl import LayerSet
from phidl import quickplot as qp
from phidl import Path, CrossSection, Device
import phidl.path as pp
import phidl.routing as pr
import phidl.geometry as pg
import phidl
from deprecated.Bertram_Functions.MyPhidlFunctions import WaveGuideMaker

'''
Here are the functions to create, using Phidl, the Schuster resonators and the Tline GDS files.
The functions are the following:
- Tline: Creates the Tline with the bondpads and the taper sections
- SchusterResonatorSmooth: Creates the Schuster resonator with the capacitors and inductors
- CapacitorSection: Creates the capacitor section of the Schuster resonator
- InductorSection: Creates the inductor section of the Schuster resonator
- SquareEtch: Creates the square etch around the Schuster resonator
- ChipResonatorsTline: Creates the chip with the Tline and the Schuster resonators
- ChipTline: Creates the chip with only the Tline

'''


def Tline(FeedlineWidth, FeedlineLength, FeedlineGap, 
          TaperLength, BondpadWidth, BondpadLength, BondpadGap,
          MWO_sim_single_res = False):
    '''
    Creates the Tline with the bondpads and the taper sections.
    Parameters:
    - FeedlineWidth: Width of the feedline
    - FeedlineLength: Length of the feedline
    - FeedlineGap: Gap of the feedline
    - TaperLength: Length of the taper sections which connect the main feedline to the bondpads
    - BondpadWidth: Width of the bondpads
    - BondpadLength: Length of the bondpads
    - BondpadGap: Gap of the bondpads
    - MWO_sim_single_res: If True, returns only the feedline without taper and bondpad.
    '''
    
    D = Device()

    #Create the feedline and bonds. First polygon is metal, second is gap.
    # Origin defined in the left edge of the feedline

    feedline, _ = WaveGuideMaker(FeedlineWidth, FeedlineGap, FeedlineLength, ['In', 'Out'])
    bondpad, _ = WaveGuideMaker(BondpadWidth, BondpadGap, BondpadLength, ['In', 'Out'])

    #References in the device and translation of bondpads
    fl = D << feedline
    bp_l = D << bondpad
    bp_l.movex(- BondpadLength - TaperLength)
    bp_r = D << bondpad
    bp_r.movex(FeedlineLength + TaperLength)

    # Taper sections
    taper_l= pr.route_smooth(fl.ports['In'], bp_l.ports['Out'], width = np.array([FeedlineWidth, BondpadWidth]), layer = 0)
    taper_l_gap = pr.route_smooth(fl.ports['In'], bp_l.ports['Out'], width = np.array([FeedlineWidth + 2*FeedlineGap, BondpadWidth + 2*BondpadGap]), layer = 1)
    taper_r = pr.route_smooth(fl.ports['Out'], bp_r.ports['In'], width = np.array([FeedlineWidth, BondpadWidth]), layer = 0)
    taper_r_gap = pr.route_smooth(fl.ports['Out'], bp_r.ports['In'], width = np.array([FeedlineWidth + 2*FeedlineGap, BondpadWidth + 2*BondpadGap]), layer = 1)
    
    #Add gap at the end of the bondpads
    width = 50
    gap_l = pg.rectangle(size = (width, BondpadWidth + 2*BondpadGap), layer = 1)
    gap_l.move((-BondpadLength - TaperLength- width, -BondpadWidth/2- BondpadGap))
    gap_r = pg.rectangle(size = (width, BondpadWidth+ 2*BondpadGap), layer = 1)
    gap_r.move((BondpadLength + TaperLength + FeedlineLength, -BondpadWidth/2 - BondpadGap))

    # We now create two Devices, one for the metallic parts 
    # and one for the gap sections (Etched parts)
    
    # Add metallic parts to a device
    Dmetal = Device('Tline_metal')
    if MWO_sim_single_res:
        Dmetal.add_polygon(fl.get_polygons()[0])
    else:
        Dmetal.add_polygon([fl.get_polygons()[0], bp_l.get_polygons()[0], 
                            bp_r.get_polygons()[0], taper_l.get_polygons(),
                            taper_r.get_polygons()])
    Dmetal = pg.union(Dmetal, layer = 0)
    
    # Add gap sections to the device
    Dgap = Device('Tline_gap')
    if MWO_sim_single_res:
        Dgap.add_polygon(fl.get_polygons()[1])
        #Add gap at the end of the bondpads
        gap_l = pg.rectangle(size = (width, FeedlineWidth + 2*FeedlineGap), layer = 1)
        gap_l.move((- width, -FeedlineWidth/2- FeedlineGap))
        gap_r = pg.rectangle(size = (width, FeedlineWidth+ 2*FeedlineGap), layer = 1)
        gap_r.move(( FeedlineLength, -FeedlineWidth/2 - FeedlineGap))
        Dgap.add_polygon([gap_l.get_polygons(), gap_r.get_polygons()])
    else:
        Dgap.add_polygon([fl.get_polygons()[1], bp_l.get_polygons()[1], 
                            bp_r.get_polygons()[1], taper_l_gap.get_polygons(),
                            taper_r_gap.get_polygons(), 
                            gap_l.get_polygons(), gap_r.get_polygons()])
    Dgap = pg.union(Dgap, layer = 1)


    return Dmetal, Dgap


def SchusterResonatorSmooth(CapacitorHorizontalLength, 
                            CapacitorVerticalLength, 
                            CapacitorWidth,
                            NumberOfBends, 
                            InductorVerticalLength,
                            InductorHorizontalLength,
                            InductorWidth,
                            InductorEndLength,
                            TaperWidth,
                            TaperLength,
                            SpacingC0, 
                            SpacingCc,
                            calib = False, 
                            cap_sim = False):
    '''
    Creates the Schuster resonator with the capacitors and inductors.
    Parameters:
    - CapacitorHorizontalLength: Horizontal length of the capacitor, which defines Cc
    - CapacitorVerticalLength: Vertical length of the capacitor, which defines C0
    - CapacitorWidth: Width of the capacitor
    - NumberOfBends: Number of bends in the inductor
    - InductorVerticalLength: Vertical length of the inductor bend
    - InductorHorizontalLength: Horizontal length of the inductor bend
    - InductorWidth: Width of the inductor
    - InductorEndLength: Length of the straight section at the end of the inductor, which connects to ground
    - TaperWidth: Width of the taper section between the inductor and the horizontal section of the capacitor
    - TaperLength: Length of the taper section between the inductor and and the horizontal section of the capacitor
    - SpacingC0: Spacing between the vertical section of the capacitor to ground. Defines C0
    - SpacingCc: Spacing between the horizontal section of the capacitor and the end of the etch box around the resonator. Defines Cc
    - calib : If True, it returns the different polygons in different layers.
    - cap_sim: If True, it returns the inductor section with a little cut to ground, so capacitive simulations can be done.
    '''
    
    D = Device()

    #Origin defined in the middle of the capacitor
    #Capacitor section

    Cap_Poly = CapacitorSection(CapacitorHorizontalLength,
                                CapacitorVerticalLength,
                                CapacitorWidth,
                                TaperWidth)
    Dcap = D << Cap_Poly
    #Inductor section
    Ind_Poly, StraightLengthInductor = InductorSection(NumberOfBends,
                                                        InductorVerticalLength,
                                                        InductorHorizontalLength,
                                                        InductorEndLength,
                                                        InductorWidth)
    Ind_Poly.move(destination = (0, -CapacitorWidth/2 - TaperLength))
    if cap_sim:
        Ind_Poly.move(destination = (0, 0.1))
    Dind = D << Ind_Poly

    # Taper section between inductor and capacitor
    route = pr.route_smooth(Dcap.ports['Midpoint'], Dind.ports['Top'], width = np.array([TaperWidth, InductorWidth]), layer = 0)
    Droute = D << route
    #Create etch box around the resonator
    Etch, ysize, xsize = SquareEtch(SpacingC0,
                                            SpacingCc,
                                            StraightLengthInductor,
                                            CapacitorHorizontalLength,
                                            CapacitorVerticalLength,
                                            CapacitorWidth, 
                                            TaperLength)

    Etch.move(origin = (0,0), destination = (-xsize/2, (-ysize + (1/2)*CapacitorWidth) + SpacingCc))
    if calib:
        return D, Etch, Dcap, Dind, Droute
    else:
        return D, Etch


def CapacitorSection(CapacitorHorizontalLength, 
                     CapacitorVerticalLength, 
                     CapacitorWidth,
                     TaperWidth):
    '''
    Creates the capacitor section of the Schuster resonator. It has a Pi shape.
    Origin defined in the middle of the horizontal length
    Parameters:
    - CapacitorHorizontalLength: Horizontal length of the capacitor, which defines Cc
    - CapacitorVerticalLength: Vertical length of the capacitor, which defines C0
    - CapacitorWidth: Width of the capacitor
    - TaperWidth: Width of the taper section between the inductor and the horizontal section of the capacitor
    '''
    #Create a path for the capacitor
    #Origin defined in the middle of the horizontal length
    points = np.array([(-CapacitorHorizontalLength/2, -CapacitorVerticalLength + CapacitorWidth/2), # CapacitorWidth/2 is the smoother radius
                       (-CapacitorHorizontalLength/2, 0),
                       (CapacitorHorizontalLength/2, 0),
                      (CapacitorHorizontalLength/2, -CapacitorVerticalLength + CapacitorWidth/2)])
    PathCapacitor = pp.smooth(points, radius = CapacitorWidth*0.5)
    PolyCapacitor = PathCapacitor.extrude(width = CapacitorWidth, simplify=1e-6, layer = 0)
    
    #Make the corners of the capacitor smooth
    SmoothLeft = pg.circle(radius = CapacitorWidth/2).move(destination = (-CapacitorHorizontalLength/2, -CapacitorVerticalLength + CapacitorWidth/2))
    SmoothRigth = pg.circle(radius = CapacitorWidth/2).move(destination = (CapacitorHorizontalLength/2, -CapacitorVerticalLength + CapacitorWidth/2))
    Smoothers = pg.boolean(SmoothLeft, SmoothRigth, operation = 'or')
    PolyCapacitor = pg.boolean(PolyCapacitor, Smoothers, operation = 'or')

    #Add port where the inductor will be connected
    PolyCapacitor.add_port(name = 'Midpoint', midpoint=(0, -CapacitorWidth/2), width = TaperWidth, orientation = 270)
    return PolyCapacitor

def InductorSection(NumberOfBends, 
                    InductorVerticalLength,
                    InductorHorizontalLength,
                    EndLength,
                    InductorWidth):
    '''
    Creates the inductor section of the Schuster resonator. It has a S-shape.
    Origin defined on the top of the inductor
    Parameters:
    - NumberOfBends: Number of bends in the inductor
    - InductorVerticalLength: Vertical length of the inductor bend. It is the vertical length of the S-shape
    - InductorHorizontalLength: Horizontal length of the inductor bend. It is the horizontal length of the S-shape
    - EndLength: Length of the straight section at the end of the inductor, which connects to ground
    - InductorWidth: Width of the inductor
    '''
    
    #Create a path for the inductor. Origin defined on the top of the inductor
    #Which will be connected to the port of the CapacitorSection
    points = [(0,0), (0, -InductorVerticalLength)]
    # Add points of the bend. Each bend is composed of 4 points that make an S-shape
    for i in range(0,NumberOfBends+1, 2):
        points.append((-InductorHorizontalLength/2, -InductorVerticalLength*(i+1)))
        points.append((-InductorHorizontalLength/2, -InductorVerticalLength*(i+2)))
        points.append((InductorHorizontalLength/2, -InductorVerticalLength*(i+2)))
        points.append((InductorHorizontalLength/2, -InductorVerticalLength*(i+3)))
        points.append((0, -InductorVerticalLength*(i+3)))
    #Add the last straight section
    if NumberOfBends==-1:
        points.append((0, -InductorVerticalLength - EndLength))
        StraightLength = InductorVerticalLength + EndLength
        PathInductor = pp.straight(length = StraightLength)
        PathInductor.rotate(-90)
    else:
        points.append((0, -InductorVerticalLength*(NumberOfBends+4) - EndLength)) #Should be +3!!!! Change later
        StraightLength = InductorVerticalLength*(NumberOfBends+4) + EndLength
        # TotalLength =  StraigthLength + InductorHorizontalLength*NumberOfBends
        PathInductor = pp.smooth(points, radius = InductorWidth*1)

    PolyInductor = PathInductor.extrude(width = InductorWidth, simplify=1e-3, layer = 0)
    PolyInductor.add_port(name = 'Top', midpoint=(0, 0), orientation = 90, width = InductorWidth)
    print(type(PolyInductor.get_polygons()[0]))
    return PolyInductor, StraightLength
    
def SquareEtch(SpacingC0, 
                SpacingCc, 
                StraightLengthInductor,
                CapacitorHorizontalLength, 
                CapacitorVerticalLength,
                CapacitorWidth,
                TaperLength):
    '''
    Creates the square etch around the Schuster resonator. The origin is defined by the capacitor section.
    Parameters:
    - SpacingC0: Spacing between the vertical section of the capacitor to ground. Defines C0
    - SpacingCc: Spacing between the horizontal section of the capacitor and the end of the etch box around the resonator. Defines Cc
    - StraightLengthInductor: Length of the straight section at the end of the inductor, which connects to ground
    - CapacitorHorizontalLength: Horizontal length of the capacitor, which defines Cc
    - CapacitorVerticalLength: Vertical length of the capacitor, which defines C0
    - CapacitorWidth: Width of the capacitor
    - TaperLength: Length of the taper section between the inductor and and the horizontal section of the capacitor
    '''
    #The origin is defined by the capacitor section

    #Square which defines the etch
    xsize = CapacitorHorizontalLength + CapacitorWidth + 2*SpacingC0
    ysize = SpacingCc 
    if StraightLengthInductor + TaperLength>CapacitorVerticalLength:
        ysize += (CapacitorWidth + StraightLengthInductor+ TaperLength) 
    else:
        raise ValueError('The vertical length of the etch box is smaller than the vertical length of the inductor')

    SquareEtch = pg.rectangle(size = (xsize, ysize), layer = 0)
    return SquareEtch, ysize, xsize


def set_layers():
    ls = LayerSet()
    ls.add_layer('Metal', gds_layer=0, color = 'red')
    ls.add_layer('EtchingBox', gds_layer=1, color = 'orange')
    ls.add_layer('Marker', gds_layer=2, color = 'green')
    ls.add_layer('NegativePlane', gds_layer=3, color = 'black')
    ls.add_layer('Resonators', gds_layer=4, color = 'blue')
    ls.add_layer('Tline', gds_layer=5, color = 'yellow')
    return ls

def chip_definition(cap_sim, Chipsize, BondpadGap, ls):
    Chip = Device('Chip')
    if cap_sim:
        FinalSpacingBondpads = 5
        Chipsize[0] = Chipsize[0] + 2*BondpadGap+50
    Chip.add_polygon(pg.rectangle(size = Chipsize,layer = ls['Metal']).get_polygons())
    Chip.move(destination = (-Chipsize[0]/2, -Chipsize[1]/2))
    return Chip

def Tline_position(FeedlineWidth, FeedlineLength, FeedlineGap, 
                    FeedlineTaperLength, BondpadWidth, BondpadLength, 
                    BondpadGap, MWO_sim_single_res, ypos,
                    cap_sim, FinalSpacingBondpads, ls):
     
    TlineMetal, TlineGap = Tline(FeedlineWidth, FeedlineLength, FeedlineGap, 
                                 FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap, MWO_sim_single_res)

    TlineMetal.movex(-FeedlineLength/2)
    TlineGap.movex(-FeedlineLength/2)
    TlineMetal.movey(ypos)
    TlineGap.movey(ypos)

    if not cap_sim:
        BondpadSpacingLeft = pg.rectangle(size = (FinalSpacingBondpads, 2*BondpadGap + BondpadWidth), layer = ls['Metal'])
        BondpadSpacingLeft.move(destination = (-FeedlineLength/2 - FeedlineTaperLength - BondpadLength - FinalSpacingBondpads, -BondpadGap - BondpadWidth/2 +ypos))
        BondpadSpacingRight = BondpadSpacingLeft.copy('BondpadSpacingRight', translation=[(FeedlineLength + 2*FeedlineTaperLength + 2*BondpadLength + FinalSpacingBondpads),0])
    else:
        BondpadSpacingRight = None
        BondpadSpacingLeft = None
    return TlineMetal, TlineGap, BondpadSpacingLeft, BondpadSpacingRight

def place_resonators(NumberOfResonators,FeedlineLength, FeedlineWidth, FeedlineGap, SeparationTlineResonator,
                     CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                    NumberOfBends, InductorVerticalLength, InductorHorizontalLength,
                    InductorWidth, InductorEndLength, TaperWidth, TaperLength,
                    SpacingC0, SpacingCc, cap_sim,
                    D_resonators, D_gap, ypos_Tline, EdgeResDistanceFactor = 0.3):
    if NumberOfResonators == 1:
        xpos = [0]
    else:
        xpos = np.linspace(-EdgeResDistanceFactor*FeedlineLength, EdgeResDistanceFactor*FeedlineLength , NumberOfResonators)
    # ypos_abs = FeedlineWidth/2 + FeedlineGap + SeparationTlineResonator + SpacingCc + CapacitorWidth/2
    sign = 1

    for i in range(NumberOfResonators):
        Resonator, Etch = SchusterResonatorSmooth(CapacitorHorizontalLength[i], 
                                                  CapacitorVerticalLength[i], 
                                                  CapacitorWidth[i],
                                                  NumberOfBends[i], 
                                                  InductorVerticalLength[i],
                                                  InductorHorizontalLength[i],
                                                  InductorWidth[i],
                                                  InductorEndLength[i],
                                                  TaperWidth[i],
                                                  TaperLength[i],
                                                  SpacingC0[i], 
                                                  SpacingCc[i], 
                                                  cap_sim = cap_sim)
        
        if sign == 1:
            Resonator.rotate(180)
            Etch.rotate(180)

        ypos = sign*(FeedlineWidth/2 + FeedlineGap + SeparationTlineResonator[i] 
                     + SpacingCc[i] + CapacitorWidth[i]/2) + ypos_Tline
        Resonator.movex(xpos[i])
        Resonator.movey(ypos)
        Etch.movex(xpos[i])
        Etch.movey(ypos)
        sign = -sign
        # D_metal.add_polygon(Resonator.get_polygons())
        # pols = Resonator.get_polygons()
        # pols.merge()
        # one_pol = pg.kl_boolean(A = Resonator.get_polygons(),B = Resonator.get_polygons(),  operation='A+B')
        # print(len(one_pol.get_polygons()))
        temp = Device('temp')
        temp.add_polygon(Resonator.get_polygons())
        # D_resonators.add_polygon(Resonator.get_polygons())
        # print(len(D_resonators.get_polygons()))
        one_pol = pg.boolean(A = temp,B = temp,  operation='A+B')
        D_resonators.add_polygon(one_pol.get_polygons())
        D_gap.add_ref(Etch)
    return D_resonators, D_gap


def FinalChipStructure(Chip, D_gap, D_metal, D_resonators, FinalSpacingBondpads, Chipsize, ls,
                       MWO_simulation, cap_sim):
    
    Ground_Plane = pg.boolean(Chip, D_gap, operation = 'not')
    EtchingBoxNegative = pg.rectangle(size = (Chipsize[0]+2*FinalSpacingBondpads, Chipsize[1] + 2*FinalSpacingBondpads), layer = ls['Metal'])
    EtchingBoxNegative.move(destination = (-Chipsize[0]/2 -FinalSpacingBondpads , -Chipsize[1]/2 - FinalSpacingBondpads)) #Center The chip
    EtchingBox = pg.boolean(EtchingBoxNegative, Chip, operation = 'not')
    Marker = pg.rectangle(size = (1000, 100), layer = ls['Marker'])
    Marker.move(destination = (-Chipsize[0]/2 + 1000, -Chipsize[1]/2 + 200))

    FinalChip = Device('FinalChip')
    FinalChip.add_polygon(Ground_Plane.get_polygons(), layer = ls['Metal'])
    FinalChip.add_polygon(D_metal.get_polygons(), layer = ls['Tline'])
    if D_resonators is not None:
        FinalChip.add_polygon(D_resonators.get_polygons(), layer = ls['Resonators'])
    if MWO_simulation or cap_sim:
        return Ground_Plane, D_metal, FinalChip
    else:
        FinalChip.add_polygon(EtchingBox.get_polygons(), layer = ls['EtchingBox'])
        FinalChip.add_polygon(Marker.get_polygons(), layer = ls['Marker'])
        
        FinalChip.add_polygon(Chip.get_polygons(), layer = ls['NegativePlane'])
    return Ground_Plane, D_metal, FinalChip


def ChipResonatorsTline(Chipsize, NumberOfResonators, SeparationTlineResonator,
                        FeedlineWidth, FeedlineLength, FeedlineGap, 
                        FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap,
                        CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                        NumberOfBends, InductorVerticalLength, InductorHorizontalLength, InductorWidth,InductorEndLength,
                        TaperWidth, TaperLength, SpacingC0, SpacingCc,
                        FinalSpacingBondpads = 100,
                        MWO_simulation = False,
                        cap_sim = False, EdgeResDistanceFactor = 0.3):
    '''
    Creates the chip with the Tline and Schuster resonators. The origin is defined in the center of the chip.
    The resonators are placed alternating in the top and bottom of the feedline.
    Parameters:
    - Chipsize: Size of the chip
    - NumberOfResonators: Number of resonators
    - SeparationTlineResonator: Separation between the feedline and the etching box of the resonators. This spacing is metal.
    - FeedlineWidth: Width of the feedline
    - FeedlineLength: Length of the feedline
    - FeedlineGap: Gap of the feedline
    - FeedlineTaperLength: Length of the taper sections which connect the main feedline to the bondpads
    - BondpadWidth: Width of the bondpads
    - BondpadLength: Length of the bondpads
    - BondpadGap: Gap of the bondpads
    - CapacitorHorizontalLength: Horizontal length of the capacitor, which defines Cc
    - CapacitorVerticalLength: Vertical length of the capacitor, which defines C0
    - CapacitorWidth: Width of the capacitor
    - NumberOfBends: Number of bends in the inductor
    - InductorVerticalLength: Vertical length of the inductor bend
    - InductorHorizontalLength: Horizontal length of the inductor bend
    - InductorWidth: Width of the inductor
    - InductorEndLength: Length of the straight section at the end of the inductor, which connects to ground
    - TaperWidth: Width of the taper section between the inductor and the horizontal section of the capacitor
    - TaperLength: Length of the taper section between the inductor and and the horizontal section of the capacitor
    - SpacingC0: Spacing between the vertical section of the capacitor to ground. Defines C0
    - SpacingCc: Spacing between the horizontal section of the capacitor and the end of the etch box around the resonator. Defines Cc
    - FinalSpacingBondpads: Spacing between the bondpads and the edge of the chip. Neeeded for dicing.
    - FluxHoles: If True, adds holes in the ground plane to avoid magnetic flux losses.
    - MWO_simulation: If True, only the metal parts are returned. If False, the metal, etching box, marker and negative plane are returned.
    - cap_sim: If True, it returns the inductor section with a little cut to ground, so capacitive simulations can be done.
    '''

    # Layers
    ls = set_layers()

    #Chip
    Chip = chip_definition(cap_sim, Chipsize, BondpadGap, ls)
    
    # Devices for the resonators and Tline metal and gap parts
    D_metal = Device('Metal')
    D_gap = Device('Gap')
    D_resonators = Device('Resonators')

    #Tline with resonators
    TlineMetal, TlineGap, BondpadSpacingLeft, BondpadSpacingRight = Tline_position(FeedlineWidth, FeedlineLength, FeedlineGap, 
                    FeedlineTaperLength, BondpadWidth, BondpadLength, 
                    BondpadGap, MWO_simulation, 0,
                    cap_sim, FinalSpacingBondpads, ls)
    
    D_metal.add_ref(TlineMetal)
    D_gap.add_ref(TlineGap)
    if  BondpadSpacingRight is not None:
        D_gap.add_polygon(BondpadSpacingLeft.get_polygons())
        D_gap.add_polygon(BondpadSpacingRight.get_polygons())

    # Resonators
    D_resonators, D_gap = place_resonators(NumberOfResonators,FeedlineLength, FeedlineWidth, FeedlineGap, SeparationTlineResonator,
                     CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                    NumberOfBends, InductorVerticalLength, InductorHorizontalLength,
                    InductorWidth, InductorEndLength, TaperWidth, TaperLength,
                    SpacingC0, SpacingCc, cap_sim,
                    D_resonators, D_gap, 0 , EdgeResDistanceFactor)
    
    #Final chip structure
    Ground_Plane, D_metal, FinalChip = FinalChipStructure(Chip, D_gap, D_metal, D_resonators, FinalSpacingBondpads, Chipsize, ls,
                       MWO_simulation, cap_sim)

    return Ground_Plane, D_metal, FinalChip

def ChipResonatorsTline_oldVersion(Chipsize, NumberOfResonators, SeparationTlineResonator,
                        FeedlineWidth, FeedlineLength, FeedlineGap, 
                        FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap,
                        CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                        NumberOfBends, InductorVerticalLength, InductorHorizontalLength, InductorWidth,InductorEndLength,
                        TaperWidth, TaperLength, SpacingC0, SpacingCc,
                        FinalSpacingBondpads = 100,
                        MWO_simulation = False,
                        cap_sim = False):
    '''
    Creates the chip with the Tline and Schuster resonators. The origin is defined in the center of the chip.
    The resonators are placed alternating in the top and bottom of the feedline.
    Parameters:
    - Chipsize: Size of the chip
    - NumberOfResonators: Number of resonators
    - SeparationTlineResonator: Separation between the feedline and the etching box of the resonators. This spacing is metal.
    - FeedlineWidth: Width of the feedline
    - FeedlineLength: Length of the feedline
    - FeedlineGap: Gap of the feedline
    - FeedlineTaperLength: Length of the taper sections which connect the main feedline to the bondpads
    - BondpadWidth: Width of the bondpads
    - BondpadLength: Length of the bondpads
    - BondpadGap: Gap of the bondpads
    - CapacitorHorizontalLength: Horizontal length of the capacitor, which defines Cc
    - CapacitorVerticalLength: Vertical length of the capacitor, which defines C0
    - CapacitorWidth: Width of the capacitor
    - NumberOfBends: Number of bends in the inductor
    - InductorVerticalLength: Vertical length of the inductor bend
    - InductorHorizontalLength: Horizontal length of the inductor bend
    - InductorWidth: Width of the inductor
    - InductorEndLength: Length of the straight section at the end of the inductor, which connects to ground
    - TaperWidth: Width of the taper section between the inductor and the horizontal section of the capacitor
    - TaperLength: Length of the taper section between the inductor and and the horizontal section of the capacitor
    - SpacingC0: Spacing between the vertical section of the capacitor to ground. Defines C0
    - SpacingCc: Spacing between the horizontal section of the capacitor and the end of the etch box around the resonator. Defines Cc
    - FinalSpacingBondpads: Spacing between the bondpads and the edge of the chip. Neeeded for dicing.
    - FluxHoles: If True, adds holes in the ground plane to avoid magnetic flux losses.
    - MWO_simulation: If True, only the metal parts are returned. If False, the metal, etching box, marker and negative plane are returned.
    - cap_sim: If True, it returns the inductor section with a little cut to ground, so capacitive simulations can be done.
    '''
    # Layers
    ls = LayerSet()
    ls.add_layer('Metal', gds_layer=0, color = 'red')
    ls.add_layer('EtchingBox', gds_layer=1, color = 'orange')
    ls.add_layer('Marker', gds_layer=2, color = 'green')
    ls.add_layer('NegativePlane', gds_layer=3, color = 'black')
    ls.add_layer('Resonators', gds_layer=4, color = 'blue')
    # if Chipsize[0] != (FeedlineLength + 2*FeedlineTaperLength + 2*BondpadLength + 2*FinalSpacingBondpads):
    #     raise ValueError(f'The chip size length ({Chipsize[0]}um) is not equal to the sum of the feedline, bondpads and spacings ({FeedlineLength + 2*FeedlineTaperLength + 2*BondpadLength + 2*FinalSpacingBondpads}um)')
    
    #Origin will be defined in the center of the chip
    Chip = Device('Chip')
    if cap_sim:
        FinalSpacingBondpads = 5
        Chipsize[0] = Chipsize[0] + 2*BondpadGap+50
    Chip.add_polygon(pg.rectangle(size = Chipsize,layer = ls['Metal']).get_polygons())
    Chip.move(destination = (-Chipsize[0]/2, -Chipsize[1]/2)) #Center The chip
    
    # Devices for the resonators and Tline metal and gap parts
    D_metal = Device('Metal')
    D_gap = Device('Gap')
    D_resonators = Device('Resonators')
    
    #Tline in the center of the chip
    TlineMetal, TlineGap = Tline(FeedlineWidth, FeedlineLength, FeedlineGap, 
                                 FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap, MWO_sim_single_res=MWO_simulation)
    TlineMetal.movex(-FeedlineLength/2)
    TlineGap.movex(-FeedlineLength/2)

    D_metal.add_ref(TlineMetal)
    D_gap.add_ref(TlineGap)

    # Resonators
    if NumberOfResonators == 1:
        xpos = [0]
    else:
        xpos = np.linspace(-0.4*FeedlineLength, 0.4*FeedlineLength , NumberOfResonators)
    # ypos_abs = FeedlineWidth/2 + FeedlineGap + SeparationTlineResonator + SpacingCc + CapacitorWidth/2
    sign = 1

    for i in range(NumberOfResonators):
        Resonator, Etch = SchusterResonatorSmooth(CapacitorHorizontalLength[i], 
                                                  CapacitorVerticalLength[i], 
                                                  CapacitorWidth[i],
                                                  NumberOfBends[i], 
                                                  InductorVerticalLength[i],
                                                  InductorHorizontalLength[i],
                                                  InductorWidth[i],
                                                  InductorEndLength[i],
                                                  TaperWidth[i],
                                                  TaperLength[i],
                                                  SpacingC0[i], 
                                                  SpacingCc[i], 
                                                  cap_sim = cap_sim)
        
        if sign == 1:
            Resonator.rotate(180)
            Etch.rotate(180)

        ypos = sign*(FeedlineWidth/2 + FeedlineGap + SeparationTlineResonator[i] 
                     + SpacingCc[i] + CapacitorWidth[i]/2)
        Resonator.movex(xpos[i])
        Resonator.movey(ypos)
        Etch.movex(xpos[i])
        Etch.movey(ypos)
        sign = -sign
        # D_metal.add_polygon(Resonator.get_polygons())
        # pols = Resonator.get_polygons()
        # pols.merge()
        # one_pol = pg.kl_boolean(A = Resonator.get_polygons(),B = Resonator.get_polygons(),  operation='A+B')
        # print(len(one_pol.get_polygons()))
        temp = Device('temp')
        temp.add_polygon(Resonator.get_polygons())
        # D_resonators.add_polygon(Resonator.get_polygons())
        # print(len(D_resonators.get_polygons()))
        one_pol = pg.boolean(A = temp,B = temp,  operation='A+B')
        D_resonators.add_polygon(one_pol.get_polygons())
        # D_resonators.add_polygon(Resonator.get_polygons())
        D_gap.add_ref(Etch)
    
    #Final spacing bondpads
    if not cap_sim:
        BondpadSpacingLeft = pg.rectangle(size = (FinalSpacingBondpads, 2*BondpadGap + BondpadWidth), layer = ls['Metal'])
        BondpadSpacingLeft.move(destination = (-FeedlineLength/2 - FeedlineTaperLength - BondpadLength - FinalSpacingBondpads, -BondpadGap - BondpadWidth/2))
        BondpadSpacingRight = BondpadSpacingLeft.copy('BondpadSpacingRight', translation=[(FeedlineLength + 2*FeedlineTaperLength + 2*BondpadLength + FinalSpacingBondpads),0])
        D_gap.add_polygon(BondpadSpacingLeft.get_polygons())
        D_gap.add_polygon(BondpadSpacingRight.get_polygons())

    #Final chip structure
    Ground_Plane = pg.boolean(Chip, D_gap, operation = 'not')
    EtchingBoxNegative = pg.rectangle(size = (Chipsize[0]+2*FinalSpacingBondpads, Chipsize[1] + 2*FinalSpacingBondpads), layer = ls['Metal'])
    EtchingBoxNegative.move(destination = (-Chipsize[0]/2 -FinalSpacingBondpads , -Chipsize[1]/2 - FinalSpacingBondpads)) #Center The chip
    EtchingBox = pg.boolean(EtchingBoxNegative, Chip, operation = 'not')
    Marker = pg.rectangle(size = (200, 200), layer = ls['Marker'])
    Marker.move(destination = (-Chipsize[0]/2 + 300, -Chipsize[1]/2 + 200))
    ## Add the holes using fill_rectangle
    # pg.xor_diff(Chip, D_gap)

    FinalChip = Device('FinalChip')
    FinalChip.add_polygon(Ground_Plane.get_polygons(), layer = ls['Metal'])
    FinalChip.add_polygon(D_metal.get_polygons(), layer = ls['Metal'])
    FinalChip.add_polygon(D_resonators.get_polygons(), layer = ls['Resonators'])
    if MWO_simulation or cap_sim:
        return Ground_Plane, D_metal, FinalChip
    else:
        FinalChip.add_polygon(EtchingBox.get_polygons(), layer = ls['EtchingBox'])
        FinalChip.add_polygon(Marker.get_polygons(), layer = ls['Marker'])
        
        FinalChip.add_polygon(Chip.get_polygons(), layer = ls['NegativePlane'])
        

        return Ground_Plane, D_metal, FinalChip

def ChipResonatorsTwoTlines(Chipsize, NumberOfResonators, SeparationTlineResonator,
                        FeedlineWidth, FeedlineLength, FeedlineGap, 
                        FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap,
                        CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                        NumberOfBends, InductorVerticalLength, InductorHorizontalLength, InductorWidth,InductorEndLength,
                        TaperWidth, TaperLength, SpacingC0, SpacingCc,
                        FinalSpacingBondpads = 100,
                        MWO_simulation = False,
                        cap_sim = False,
                        ypos_tlines = [0,0],
                        EdgeResDistanceFactor = 0.3):
    
    '''
    Creates the chip with two Tlines - the first one with resonators and the second one without resonators.
    The origin is defined in the center of the chip.

    Parameters:
    (same as ChipTlineResonators)
    - ypos_tlines: Y position of the Tlines. The first element is the Tline with resonators and the second element is the Tline without resonators.
    '''
    
    # Layers
    ls = set_layers()

    #Chip
    Chip = chip_definition(cap_sim, Chipsize, BondpadGap, ls)
    
    # Devices for the resonators and Tline metal and gap parts
    D_metal = Device('Metal')
    D_gap = Device('Gap')
    D_resonators = Device('Resonators')
    
    #Tline with resonators
    TlineMetal, TlineGap, BondpadSpacingLeft, BondpadSpacingRight = Tline_position(FeedlineWidth, FeedlineLength, FeedlineGap, 
                    FeedlineTaperLength, BondpadWidth, BondpadLength, 
                    BondpadGap, MWO_simulation, ypos_tlines[0],
                    cap_sim, FinalSpacingBondpads, ls)

    D_metal.add_ref(TlineMetal)
    D_gap.add_ref(TlineGap)
    if  BondpadSpacingRight is not None:
        D_gap.add_polygon(BondpadSpacingLeft.get_polygons())
        D_gap.add_polygon(BondpadSpacingRight.get_polygons())

    #Tline without resonators
    TlineMetal_wo, TlineGap_wo, BondpadSpacingLeft_wo, BondpadSpacingRight_wo = Tline_position(FeedlineWidth, FeedlineLength, FeedlineGap, 
                    FeedlineTaperLength, BondpadWidth, BondpadLength, 
                    BondpadGap, MWO_simulation, ypos_tlines[1],
                    cap_sim, FinalSpacingBondpads, ls)

    D_metal.add_ref(TlineMetal_wo)
    D_gap.add_ref(TlineGap_wo)

    if  BondpadSpacingRight_wo is not None:
        D_gap.add_polygon(BondpadSpacingLeft_wo.get_polygons())
        D_gap.add_polygon(BondpadSpacingRight_wo.get_polygons())

    #Resonators
    D_resonators, D_gap = place_resonators(NumberOfResonators,FeedlineLength, FeedlineWidth, FeedlineGap, SeparationTlineResonator,
                     CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                    NumberOfBends, InductorVerticalLength, InductorHorizontalLength,
                    InductorWidth, InductorEndLength, TaperWidth, TaperLength,
                    SpacingC0, SpacingCc, cap_sim,
                    D_resonators, D_gap, ypos_tlines[0], EdgeResDistanceFactor)
    
    #Final chip structure
    Ground_Plane, D_metal, FinalChip = FinalChipStructure(Chip, D_gap, D_metal, D_resonators, FinalSpacingBondpads, Chipsize, ls,
                       MWO_simulation, cap_sim)

    return Ground_Plane, D_metal, FinalChip
    


def ChipResonatorsThreeTlines(Chipsize, NumberOfResonators, SeparationTlineResonator,
                        FeedlineWidth, FeedlineLength, FeedlineGap, 
                        FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap,
                        CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                        NumberOfBends, InductorVerticalLength, InductorHorizontalLength, InductorWidth,InductorEndLength,
                        TaperWidth, TaperLength, SpacingC0, SpacingCc,
                        FinalSpacingBondpads = 100,
                        MWO_simulation = False,
                        cap_sim = False,
                        ypos_tlines = [0,0,0]):
    '''
    Creates the chip with three Tlines - the first two with resonators and the third one without resonators.
    The origin is defined in the center of the chip.

    Parameters:
    (same as ChipTlineResonators)
    - ypos_tlines: Y position of the Tlines. The first and second elements are the Tline with resonators and the third element is the Tline without resonators.
    '''
    # Layers
    ls = set_layers()

    #Chip
    Chip = chip_definition(cap_sim, Chipsize, BondpadGap, ls)
    
    # Devices for the resonators and Tline metal and gap parts
    D_metal = Device('Metal')
    D_gap = Device('Gap')
    D_resonators = Device('Resonators')
    
    #Tline with resonators
    TlineMetal, TlineGap, BondpadSpacingLeft, BondpadSpacingRight = Tline_position(FeedlineWidth, FeedlineLength, FeedlineGap, 
                    FeedlineTaperLength, BondpadWidth, BondpadLength, 
                    BondpadGap, MWO_simulation, ypos_tlines[0],
                    cap_sim, FinalSpacingBondpads, ls)

    D_metal.add_ref(TlineMetal)
    D_gap.add_ref(TlineGap)
    if  BondpadSpacingRight is not None:
        D_gap.add_polygon(BondpadSpacingLeft.get_polygons())
        D_gap.add_polygon(BondpadSpacingRight.get_polygons())

    #Tline with resonators2
    TlineMetal2, TlineGap2, BondpadSpacingLeft2, BondpadSpacingRight2 = Tline_position(FeedlineWidth, FeedlineLength, FeedlineGap, 
                    FeedlineTaperLength, BondpadWidth, BondpadLength, 
                    BondpadGap, MWO_simulation, ypos_tlines[1],
                    cap_sim, FinalSpacingBondpads, ls)

    D_metal.add_ref(TlineMetal2)
    D_gap.add_ref(TlineGap2)
    if  BondpadSpacingRight2 is not None:
        D_gap.add_polygon(BondpadSpacingLeft2.get_polygons())
        D_gap.add_polygon(BondpadSpacingRight2.get_polygons())



    #Tline without resonators
    TlineMetal_wo, TlineGap_wo, BondpadSpacingLeft_wo, BondpadSpacingRight_wo = Tline_position(FeedlineWidth, FeedlineLength, FeedlineGap, 
                    FeedlineTaperLength, BondpadWidth, BondpadLength, 
                    BondpadGap, MWO_simulation, ypos_tlines[2],
                    cap_sim, FinalSpacingBondpads, ls)

    D_metal.add_ref(TlineMetal_wo)
    D_gap.add_ref(TlineGap_wo)

    if  BondpadSpacingRight_wo is not None:
        D_gap.add_polygon(BondpadSpacingLeft_wo.get_polygons())
        D_gap.add_polygon(BondpadSpacingRight_wo.get_polygons())

    #Resonators
    D_resonators, D_gap = place_resonators(NumberOfResonators,FeedlineLength, FeedlineWidth, FeedlineGap, SeparationTlineResonator,
                     CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                    NumberOfBends, InductorVerticalLength, InductorHorizontalLength,
                    InductorWidth, InductorEndLength, TaperWidth, TaperLength,
                    SpacingC0, SpacingCc, cap_sim,
                    D_resonators, D_gap, ypos_tlines[0] )
    
    #Resonators2
    D_resonators, D_gap = place_resonators(NumberOfResonators,FeedlineLength, FeedlineWidth, FeedlineGap, SeparationTlineResonator,
                     CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                    NumberOfBends, InductorVerticalLength, InductorHorizontalLength,
                    InductorWidth, InductorEndLength, TaperWidth, TaperLength,
                    SpacingC0, SpacingCc, cap_sim,
                    D_resonators, D_gap, ypos_tlines[1] )
    
    #Final chip structure
    Ground_Plane, D_metal, FinalChip = FinalChipStructure(Chip, D_gap, D_metal, D_resonators, FinalSpacingBondpads, Chipsize, ls,
                       MWO_simulation, cap_sim)

    return Ground_Plane, D_metal, FinalChip

def include_Tline(FeedlineWidth, FeedlineLength, FeedlineGap, 
                    FeedlineTaperLength, BondpadWidth, BondpadLength, 
                    BondpadGap, MWO_simulation, ypos_tlines,
                    cap_sim, FinalSpacingBondpads, ls, D_metal, D_gap, BondpadSpacingRight, BondpadSpacingLeft):
    
    TlineMetal, TlineGap, BondpadSpacingLeft, BondpadSpacingRight = Tline_position(FeedlineWidth, FeedlineLength, FeedlineGap, 
                    FeedlineTaperLength, BondpadWidth, BondpadLength, 
                    BondpadGap, MWO_simulation, ypos_tlines,
                    cap_sim, FinalSpacingBondpads, ls)

    D_metal.add_ref(TlineMetal)
    D_gap.add_ref(TlineGap)
    if  BondpadSpacingRight is not None:
        D_gap.add_polygon(BondpadSpacingLeft.get_polygons())
        D_gap.add_polygon(BondpadSpacingRight.get_polygons())
    return D_metal, D_gap

def ChipAndersen(Chipsize, NumberOfResonators, SeparationTlineResonator,
                        FeedlineWidth, FeedlineLength, FeedlineGap, 
                        FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap,
                        CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                        NumberOfBends, InductorVerticalLength, InductorHorizontalLength, InductorWidth,InductorEndLength,
                        TaperWidth, TaperLength, SpacingC0, SpacingCc,
                        FinalSpacingBondpads = 100,
                        MWO_simulation = False,
                        cap_sim = False,
                        ypos_tlines = [0,0,0]):
    ''' Creates the chip for 9x9 chip with 5 Tlines. 1 without resonators, 2 with a set of resonators and other two with another set of resonators the Andersen resonators and Tline.
      The origin is defined in the center of the chip.
    The resonators are placed alternating in the top and bottom of the feedline.
    The parameters of the resonators are tuples with the parameters of the two sets of resonators.
    '''
    # Layers
    ls = set_layers()

    #Chip
    Chip = chip_definition(cap_sim, Chipsize, BondpadGap, ls)
    
    # Devices for the resonators and Tline metal and gap parts
    D_metal = Device('Metal')
    D_gap = Device('Gap')
    D_resonators = Device('Resonators')
    
    #Tlines 
    for ypos in ypos_tlines:
        D_metal, D_gap = include_Tline(FeedlineWidth, FeedlineLength, FeedlineGap, 
                        FeedlineTaperLength, BondpadWidth, BondpadLength, 
                        BondpadGap, MWO_simulation, ypos,
                        cap_sim, FinalSpacingBondpads, ls, D_metal, D_gap, BondpadSpacingRight = None, BondpadSpacingLeft = None)
    


    # Set of resonators 1
    D_resonators, D_gap = place_resonators(NumberOfResonators,FeedlineLength, FeedlineWidth, FeedlineGap, SeparationTlineResonator[0],
                     CapacitorHorizontalLength[0], CapacitorVerticalLength[0], CapacitorWidth[0],
                    NumberOfBends[0], InductorVerticalLength[0], InductorHorizontalLength[0],
                    InductorWidth[0], InductorEndLength[0], TaperWidth[0], TaperLength[0],
                    SpacingC0[0], SpacingCc[0], cap_sim,
                    D_resonators, D_gap, ypos_tlines[0] )
    
    # Set of resonators 1 again
    D_resonators, D_gap = place_resonators(NumberOfResonators,FeedlineLength, FeedlineWidth, FeedlineGap, SeparationTlineResonator[0],
                     CapacitorHorizontalLength[0], CapacitorVerticalLength[0], CapacitorWidth[0],
                    NumberOfBends[0], InductorVerticalLength[0], InductorHorizontalLength[0],
                    InductorWidth[0], InductorEndLength[0], TaperWidth[0], TaperLength[0],
                    SpacingC0[0], SpacingCc[0], cap_sim,
                    D_resonators, D_gap, ypos_tlines[1] )
    
    # Set of resonators 2
    D_resonators, D_gap = place_resonators(NumberOfResonators,FeedlineLength, FeedlineWidth, FeedlineGap, SeparationTlineResonator[1],
                     CapacitorHorizontalLength[1], CapacitorVerticalLength[1], CapacitorWidth[1],
                    NumberOfBends[1], InductorVerticalLength[1], InductorHorizontalLength[1],
                    InductorWidth[1], InductorEndLength[1], TaperWidth[1], TaperLength[1],
                    SpacingC0[1], SpacingCc[1], cap_sim,
                    D_resonators, D_gap, ypos_tlines[2] )
    
    # Set of resonators 2 again
    D_resonators, D_gap = place_resonators(NumberOfResonators,FeedlineLength, FeedlineWidth, FeedlineGap, SeparationTlineResonator[1],
                     CapacitorHorizontalLength[1], CapacitorVerticalLength[1], CapacitorWidth[1],
                    NumberOfBends[1], InductorVerticalLength[1], InductorHorizontalLength[1],
                    InductorWidth[1], InductorEndLength[1], TaperWidth[1], TaperLength[1],
                    SpacingC0[1], SpacingCc[1], cap_sim,
                    D_resonators, D_gap, ypos_tlines[3] )

    #Final chip structure
    Ground_Plane, D_metal, FinalChip = FinalChipStructure(Chip, D_gap, D_metal, D_resonators, FinalSpacingBondpads, Chipsize, ls,
                       MWO_simulation, cap_sim)

    return Ground_Plane, D_metal, FinalChip


def ChipTline(Chipsize,
            FeedlineWidth, FeedlineLength, FeedlineGap, 
            FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap):
    '''
    Creates the chip with only the Tline. The origin is defined in the center of the chip.
    Parameters:
    - Chipsize: Size of the chip
    - FeedlineWidth: Width of the feedline
    - FeedlineLength: Length of the feedline
    - FeedlineGap: Gap of the feedline
    - FeedlineTaperLength: Length of the taper sections which connect the main feedline to the bondpads
    - BondpadWidth: Width of the bondpads
    - BondpadLength: Length of the bondpads
    - BondpadGap: Gap of the bondpads
    '''

    
    # Layers
    ls = set_layers()

    #Chip
    Chip = chip_definition(False, Chipsize, BondpadGap, ls)
    
    # Devices for the resonators and Tline metal and gap parts
    D_metal = Device('Metal')
    D_gap = Device('Gap')
    D_resonators = Device('Resonators')

    #Tline with resonators
    TlineMetal, TlineGap, BondpadSpacingLeft, BondpadSpacingRight = Tline_position(FeedlineWidth, FeedlineLength, FeedlineGap, 
                    FeedlineTaperLength, BondpadWidth, BondpadLength, 
                    BondpadGap, False, 0,
                    False, 100, ls)
    
    D_metal.add_ref(TlineMetal)
    D_gap.add_ref(TlineGap)
    if  BondpadSpacingRight is not None:
        D_gap.add_polygon(BondpadSpacingLeft.get_polygons())
        D_gap.add_polygon(BondpadSpacingRight.get_polygons())

    
    #Final chip structure
    Ground_Plane, D_metal, FinalChip = FinalChipStructure(Chip, D_gap, D_metal, None, 100, Chipsize, ls,
                       False, False)

    return Ground_Plane, D_metal, FinalChip



    # # Layers
    # ls = set_layers()
    # # ls.add_layer('Ground', gds_layer=0, color = 'red')
    # # ls.add_layer('Metal', gds_layer=1, color = 'blue')
    # # ls.add_layer('Negative plane', gds_layer=2, color = 'black')
    # #Origin will be defined in the center of the chip

    # Chip = Device('Chip')
    # Chip.add_polygon(pg.rectangle(size = Chipsize,layer = ls['Ground']).get_polygons())
    # Chip.move(destination = (-Chipsize[0]/2, -Chipsize[1]/2)) #Center The chip
    
    # # Devices for the resonators and Tline metal and gap parts
    # D_metal = Device('Metal')
    # D_gap = Device('Gap')
    
    # #Tline in the center of the chip
    # TlineMetal, TlineGap = Tline(FeedlineWidth, FeedlineLength, FeedlineGap, 
    #                              FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap)
    # TlineMetal.movex(-FeedlineLength/2)
    # TlineGap.movex(-FeedlineLength/2)

    # D_metal.add_ref(TlineMetal)
    # D_gap.add_ref(TlineGap)

    # Ground_Plane = pg.boolean(Chip, D_gap, operation = 'not')
    # FinalChip = Device('FinalChip')
    # FinalChip.add_polygon(Ground_Plane.get_polygons(), layer = ls['Ground'])
    # FinalChip.add_polygon(D_metal.get_polygons(), layer = ls['Ground'])
    # FinalChip.add_polygon(Chip.get_polygons(), layer = ls['Negative plane'])

    # return Ground_Plane, D_metal, FinalChip
