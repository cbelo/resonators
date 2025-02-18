
import ScriptEnv
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.NewProject()
oProject.InsertDesign("Q3D Extractor", "Q3DDesign1", "", "")
oDesign = oProject.SetActiveDesign("Q3DDesign1")
oDesign.RenameDesignInstance("Q3DDesign1", "scripttest")
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.ImportGDSII(
    [
        "NAME:options",
        "FileName:="        , "C:/Users/localadmin/Documents/Ansoft/julian_gds/tooth_testchip_zerovl_precision10_extraDX_extraDiaglen.gds",
        "FlattenHierarchy:="    , True,
        [
            "NAME:LayerMap",
            [
                "NAME:LayerMapInfo",
                "LayerNum:="        , 0,
                "DestLayer:="        , "Signal0",
                "layer_type:="        , "signal"
            ],
            [
                "NAME:LayerMapInfo",
                "LayerNum:="        , 4,
                "DestLayer:="        , "Signal4",
                "layer_type:="        , "signal"
            ]
        ],
        "OrderMap:="        , [            "entry:="        , [                "order:="        , 0,                "layer:="        , "Signal0"],            "entry:="        , [                "order:="        , 1,                "layer:="        , "Signal4"]]
    ])

# oDesign = oProject.SetActiveDesign("scripttest")
# oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.CreateBox(
    [
        "NAME:BoxParameters",
        "XPosition:="        , "-2.5125mm",
        "YPosition:="        , "-1.5mm",
        "ZPosition:="        , "0mm",
        "XSize:="        , "5.025mm",
        "YSize:="        , "3mm",
        "ZSize:="        , "-0.6mm"
    ], 
    [
        "NAME:Attributes",
        "Name:="        , "Box1",
        "Flags:="        , "",
        "Color:="        , "(143 175 143)",
        "Transparency:="    , 0,
        "PartCoordinateSystem:=", "Global",
        "UDMId:="        , "",
        "MaterialValue:="    , "\"copper\"",
        "SurfaceMaterialValue:=", "\"\"",
        "SolveInside:="        , False,
        "ShellElement:="    , False,
        "ShellElementThickness:=", "0mm",
        "ReferenceTemperature:=", "20cel",
        "IsMaterialEditable:="    , True,
        "IsSurfaceMaterialEditable:=", True,
        "UseMaterialAppearance:=", False,
        "IsLightweight:="    , False
    ])
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.AssignMaterial(
    [
        "NAME:Selections",
        "AllowRegionDependentPartSelectionForPMLCreation:=", True,
        "AllowRegionSelectionForPMLCreation:=", True,
        "Selections:="        , "Box1"
    ], 
    [
        "NAME:Attributes",
        "MaterialValue:="    , "\"silicon\"",
        "SolveInside:="        , False,
        "ShellElement:="    , False,
        "ShellElementThickness:=", "nan ",
        "ReferenceTemperature:=", "nan ",
        "IsMaterialEditable:="    , True,
        "IsSurfaceMaterialEditable:=", True,
        "UseMaterialAppearance:=", False,
        "IsLightweight:="    , False
    ])
oModule = oDesign.GetModule("BoundarySetup")
oModule.AssignThinConductor(
    [
        "NAME:ThinCond1",
        "Objects:="        , ["Signal0_1","Signal0_2","Signal0_3","Signal0_4","Signal4_5","Signal4_6","Signal4_7","Signal4_8","Signal4_9","Signal4_10","Signal4_11","Signal4_12","Signal4_13","Signal4_14","Signal4_15","Signal4_16","Signal4_17","Signal4_18","Signal4_19","Signal4_20","Signal4_21","Signal4_22","Signal4_23","Signal4_24","Signal4_25","Signal4_26","Signal4_27","Signal4_28"],
        "Material:="        , "Copper",
        "Thickness:="        , "100nm"
    ])
oModule.AutoIdentifyNets()
oModule = oDesign.GetModule("AnalysisSetup")
oModule.InsertSetup("Matrix", 
    [
        "NAME:Setup1",
        "AdaptiveFreq:="    , "5GHz",
        "SaveFields:="        , False,
        "Enabled:="        , True,
        [
            "NAME:Cap",
            "MaxPass:="        , 10,
            "MinPass:="        , 1,
            "MinConvPass:="        , 1,
            "PerError:="        , 1,
            "PerRefine:="        , 30,
            "AutoIncreaseAccuracy:=", True,
            "AccuracyLevel:="    , "High",
            "Solver Type:="        , "Iterative"
        ],
        "EnableTransitionRegionSolve:=", False,
        "ErrorValue:="        , "0.1"
    ])
oProject.Save()
oDesign.Analyze("Setup1")
oModule = oDesign.GetModule("BoundarySetup")
oModule.RenameBoundary("Signal4_5", "res2")
oModule.RenameBoundary("Signal0_1", "gnd")
oModule.RenameBoundary("Signal4_9", "res4")
oModule.RenameBoundary("Signal4_21", "res3")
oModule = oDesign.GetModule("ReportSetup")
oModule.CreateReport("C Matrix Table 1", "Matrix", "Data Table", "Setup1 : LastAdaptive", 
    [
        "Context:="        , "Original"
    ], 
    [
        "Freq:="        , ["All"],
        "extent_x_pos:="    , ["Nominal"],
        "extent_y_pos:="    , ["Nominal"],
        "extent_x_size:="    , ["Nominal"],
        "extent_y_size:="    , ["Nominal"],
        "Signal1_lower_elevation:=", ["Nominal"],
        "Signal2_lower_elevation:=", ["Nominal"],
        "extent_x_pos_1:="    , ["Nominal"],
        "extent_y_pos_1:="    , ["Nominal"],
        "extent_x_size_1:="    , ["Nominal"],
        "extent_y_size_1:="    , ["Nominal"],
        "Signal1_lower_elevation_1:=", ["Nominal"],
        "Signal2_lower_elevation_1:=", ["Nominal"],
        "extent_x_pos_2:="    , ["Nominal"],
        "extent_y_pos_2:="    , ["Nominal"],
        "extent_x_size_2:="    , ["Nominal"],
        "extent_y_size_2:="    , ["Nominal"],
        "Signal0_lower_elevation:=", ["Nominal"],
        "Signal4_lower_elevation:=", ["Nominal"]
    ], 
    [
        "X Component:="        , "Freq",
        "Y Component:="        , ["C(res2,gnd)","C(res3,gnd)","C(res4,gnd)","C(res3,res2)","C(res4,res2)","C(res4,res3)"]
    ])
