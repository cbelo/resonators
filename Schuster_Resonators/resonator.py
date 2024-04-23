import numpy as np
from phidl import LayerSet
from phidl import quickplot as qp
from phidl import Path, CrossSection, Device
import phidl.path as pp
import phidl.routing as pr
import phidl.geometry as pg
import phidl
from Bertram_Functions.MyPhidlFunctions import WaveGuideMaker

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
          TaperLength, BondpadWidth, BondpadLength, BondpadGap):
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
    Dmetal.add_polygon([fl.get_polygons()[0], bp_l.get_polygons()[0], 
                        bp_r.get_polygons()[0], taper_l.get_polygons(),
                          taper_r.get_polygons()])
    Dmetal = pg.union(Dmetal, layer = 0)
    
    # Add gap sections to the device
    Dgap = Device('Tline_gap')
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
                            SpacingCc):
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
    PathCapacitor = pp.smooth(points, radius = CapacitorWidth)
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
    points.append((0, -InductorVerticalLength*(NumberOfBends+4) - EndLength))
    

    StraightLength = InductorVerticalLength*(NumberOfBends+4) + EndLength
    # TotalLength =  StraigthLength + InductorHorizontalLength*NumberOfBends
    PathInductor = pp.smooth(points, radius = InductorWidth*2)
    PolyInductor = PathInductor.extrude(width = InductorWidth, simplify=1e-1, layer = 0)
    PolyInductor.add_port(name = 'Top', midpoint=(0, 0), orientation = 90, width = InductorWidth)
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
    if StraightLengthInductor>CapacitorVerticalLength:
        ysize += (CapacitorWidth + StraightLengthInductor+ TaperLength) 
    else:
        ysize += CapacitorVerticalLength
    SquareEtch = pg.rectangle(size = (xsize, ysize), layer = 0)
    return SquareEtch, ysize, xsize
                    
def ChipResonatorsTline(Chipsize, NumberOfResonators, SeparationTlineResonator,
                        FeedlineWidth, FeedlineLength, FeedlineGap, 
                        FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap,
                        CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                        NumberOfBends, InductorVerticalLength, InductorHorizontalLength, InductorWidth,InductorEndLength,
                        TaperWidth, TaperLength, SpacingC0, SpacingCc,
                        FinalSpacingBondpads = 100):
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
    '''
    # Layers
    ls = LayerSet()
    ls.add_layer('Ground', gds_layer=0, color = 'red')
    ls.add_layer('Metal', gds_layer=1, color = 'blue')
    if Chipsize[0] != (FeedlineLength + 2*FeedlineTaperLength + 2*BondpadLength + 2*FinalSpacingBondpads):
        raise ValueError(f'The chip size length ({Chipsize[0]}um) is not equal to the sum of the feedline, bondpads and spacings ({FeedlineLength + 2*FeedlineTaperLength + 2*BondpadLength + 2*FinalSpacingBondpads}um)')
    
    #Origin will be defined in the center of the chip
    Chip = Device('Chip')
    Chip.add_polygon(pg.rectangle(size = Chipsize,layer = ls['Ground']).get_polygons())
    Chip.move(destination = (-Chipsize[0]/2, -Chipsize[1]/2)) #Center The chip
    
    # Devices for the resonators and Tline metal and gap parts
    D_metal = Device('Metal')
    D_gap = Device('Gap')
    
    #Tline in the center of the chip
    TlineMetal, TlineGap = Tline(FeedlineWidth, FeedlineLength, FeedlineGap, 
                                 FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap)
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
                                                  SpacingCc[i])
        
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
        D_metal.add_polygon(Resonator.get_polygons())
        D_gap.add_ref(Etch)
    
    #Final spacing bondpads
    BondpadSpacingLeft = pg.rectangle(size = (FinalSpacingBondpads, 2*BondpadGap + BondpadWidth), layer = ls['Ground'])
    BondpadSpacingLeft.move(destination = (-FeedlineLength/2 - FeedlineTaperLength - BondpadLength - FinalSpacingBondpads, -BondpadGap - BondpadWidth/2))
    BondpadSpacingRight = BondpadSpacingLeft.copy('BondpadSpacingRight', translation=[(FeedlineLength + 2*FeedlineTaperLength + 2*BondpadLength + FinalSpacingBondpads),0])
    D_gap.add_polygon(BondpadSpacingLeft.get_polygons())
    D_gap.add_polygon(BondpadSpacingRight.get_polygons())

    #Final chip structure
    Ground_Plane = pg.boolean(Chip, D_gap, operation = 'not')

    ## Add the holes using fill_rectangle
    # pg.xor_diff(Chip, D_gap)

    FinalChip = Device('FinalChip')
    FinalChip.add_polygon(Ground_Plane.get_polygons(), layer = ls['Ground'])
    FinalChip.add_polygon(D_metal.get_polygons(), layer = ls['Metal'])

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
    ls = LayerSet()
    ls.add_layer('Ground', gds_layer=0, color = 'red')
    ls.add_layer('Metal', gds_layer=1, color = 'blue')

    #Origin will be defined in the center of the chip
    Chip = Device('Chip')
    Chip.add_polygon(pg.rectangle(size = Chipsize,layer = ls['Ground']).get_polygons())
    Chip.move(destination = (-Chipsize[0]/2, -Chipsize[1]/2)) #Center The chip
    
    # Devices for the resonators and Tline metal and gap parts
    D_metal = Device('Metal')
    D_gap = Device('Gap')
    
    #Tline in the center of the chip
    TlineMetal, TlineGap = Tline(FeedlineWidth, FeedlineLength, FeedlineGap, 
                                 FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap)
    TlineMetal.movex(-FeedlineLength/2)
    TlineGap.movex(-FeedlineLength/2)

    D_metal.add_ref(TlineMetal)
    D_gap.add_ref(TlineGap)

    
    Ground_Plane = pg.boolean(Chip, D_gap, operation = 'not')
    FinalChip = Device('FinalChip')
    FinalChip.add_polygon(Ground_Plane.get_polygons(), layer = ls['Ground'])
    FinalChip.add_polygon(D_metal.get_polygons(), layer = ls['Metal'])

    return Ground_Plane, D_metal, FinalChip
