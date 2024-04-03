import numpy as np
from phidl import LayerSet
from phidl import quickplot as qp
from phidl import Path, CrossSection, Device
import phidl.path as pp
import phidl.routing as pr
import phidl.geometry as pg
import phidl
from MyPhidlFunctions import WaveGuideMaker


def Tline(FeedlineWidth, FeedlineLength, FeedlineGap, 
          TaperLength, BondpadWidth, BondpadLength, BondpadGap):
    
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
                            TaperWidth,
                            TaperLength,
                            SpacingC0, 
                            SpacingCc):
    
    D = Device()

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
                                                        CapacitorWidth,
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
    
    Etch.move(destination = (-xsize/2, -ysize + (3/2)*CapacitorWidth))
    
    return D, Etch


def CapacitorSection(CapacitorHorizontalLength, 
                     CapacitorVerticalLength, 
                     CapacitorWidth,
                     TaperWidth):
    
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
    #The prigin is defined by the capacitor section

    #Square which defines the etch
    xsize = CapacitorHorizontalLength + CapacitorWidth + 2*SpacingC0
    ysize = SpacingCc + TaperLength
    if StraightLengthInductor>CapacitorVerticalLength:
        ysize += (CapacitorWidth + StraightLengthInductor) 
    else:
        ysize += CapacitorVerticalLength
    SquareEtch = pg.rectangle(size = (xsize, ysize), layer = 0)
    return SquareEtch, ysize, xsize
                    
def ChipResonatorsTline(Chipsize, NumberOfResonators, SeparationTlineResonator,
                        FeedlineWidth, FeedlineLength, FeedlineGap, 
                        FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap,
                        CapacitorHorizontalLength, CapacitorVerticalLength, CapacitorWidth,
                        NumberOfBends, InductorVerticalLength, InductorHorizontalLength, InductorWidth,
                        TaperWidth, TaperLength, SpacingC0, SpacingCc):
    
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

    # Resonators
    xpos = np.linspace(-FeedlineLength/2 + 1.2*CapacitorHorizontalLength[0], FeedlineLength/2 - 1.2*CapacitorHorizontalLength[-1] , NumberOfResonators)
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
                                                  TaperWidth[i],
                                                  TaperLength[i],
                                                  SpacingC0[i], 
                                                  SpacingCc[i])
        
        if sign == 1:
            Resonator.rotate(180)
            Etch.rotate(180)

        ypos = sign*(FeedlineWidth/2 + FeedlineGap + SeparationTlineResonator 
                     + SpacingCc[i] + CapacitorWidth[i]/2)
        Resonator.movex(xpos[i])
        Resonator.movey(ypos)
        Etch.movex(xpos[i])
        Etch.movey(ypos)
        sign = -sign
        D_metal.add_polygon(Resonator.get_polygons())
        D_gap.add_ref(Etch)
    
    Ground_Plane = pg.boolean(Chip, D_gap, operation = 'not')
    # pg.xor_diff(Chip, D_gap)

    FinalChip = Device('FinalChip')
    FinalChip.add_polygon(Ground_Plane.get_polygons(), layer = ls['Ground'])
    FinalChip.add_polygon(D_metal.get_polygons(), layer = ls['Metal'])

    return Ground_Plane, D_metal, FinalChip


def ChipTline(Chipsize,
            FeedlineWidth, FeedlineLength, FeedlineGap, 
            FeedlineTaperLength, BondpadWidth, BondpadLength, BondpadGap,
):
    
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
