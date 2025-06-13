# ----------------------------------------------
# Script Recorded by Ansys Electronics Desktop Version 2024.2.0
# 17:01:29  Feb 18, 2025
# ----------------------------------------------
import ScriptEnv
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.SetActiveProject("test1")
oDesign = oProject.SetActiveDesign("Q3DDesign1")
oEditor = oDesign.SetActiveEditor("3D Modeler")
ysize = 1.2
ypos = -ysize/2
chipsize_in = 1100
BondpadGap = 10.5747
xsize = chipsize_in*0.001 + 2*BondpadGap*0.001 + 0.05
xpos = -xsize/2
xpos = '{:.6f}mm'.format(xpos)
ypos = '{:.6f}mm'.format(ypos)
xsize = '{:.6f}mm'.format(xsize)
ysize = '{:.6f}mm'.format(ysize)
oEditor.CreateBox(
	[
		"NAME:BoxParameters",
		"XPosition:="		, xpos,
		"YPosition:="		, ypos,
		"ZPosition:="		, "0mm",
		"XSize:="		, xsize,
		"YSize:="		, ysize,
		"ZSize:="		, "-0.1mm"
	], 
	[
		"NAME:Attributes",
		"Name:="		, "Box1",
		"Flags:="		, "",
		"Color:="		, "(143 175 143)",
		"Transparency:="	, 0,
		"PartCoordinateSystem:=", "Global",
		"UDMId:="		, "",
		"MaterialValue:="	, "\"copper\"",
		"SurfaceMaterialValue:=", "\"\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "0mm",
		"ReferenceTemperature:=", "20cel",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])
oEditor.AssignMaterial(
	[
		"NAME:Selections",
		"AllowRegionDependentPartSelectionForPMLCreation:=", True,
		"AllowRegionSelectionForPMLCreation:=", True,
		"Selections:="		, "Box1"
	], 
	[
		"NAME:Attributes",
		"MaterialValue:="	, "\"SiGe\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "nan ",
		"ReferenceTemperature:=", "nan ",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])

oEditor.ThickenSheet(
	[
		"NAME:Selections",
		"Selections:="		, "Signal0_1,Signal0_2,Signal4_3",
		"NewPartsModelFlag:="	, "Model"
	], 
	[
		"NAME:SheetThickenParameters",
		"Thickness:="		, "-0.05mm",
		"BothSides:="		, False,
		[
			"NAME:ThickenAdditionalInfo",
			[
				"NAME:ShellThickenDirectionInfo",
				"SampleFaceID:="	, 18,
				"ComponentSense:="	, True,
				[
					"NAME:PointOnSampleFace",
					"X:="			, "-0.320404859475mm",
					"Y:="			, "-0.007136581mm",
					"Z:="			, "0mm"
				],
				[
					"NAME:DirectionAtPoint",
					"X:="			, "0mm",
					"Y:="			, "0mm",
					"Z:="			, "1mm"
				]
			]
		]
	])
oEditor.AssignMaterial(
	[
		"NAME:Selections",
		"AllowRegionDependentPartSelectionForPMLCreation:=", True,
		"AllowRegionSelectionForPMLCreation:=", True,
		"Selections:="		, "Signal0_1,Signal0_2,Signal4_3"
	], 
	[
		"NAME:Attributes",
		"MaterialValue:="	, "\"platinum\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "nan ",
		"ReferenceTemperature:=", "nan ",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])
oEditor.Subtract(
	[
		"NAME:Selections",
		"Blank Parts:="		, "Box1",
		"Tool Parts:="		, "Signal0_1,Signal0_2,Signal4_3"
	], 
	[
		"NAME:SubtractParameters",
		"KeepOriginals:="	, True,
		"TurnOnNBodyBoolean:="	, True
	])
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.ChangeProperty(
	[
		"NAME:AllTabs",
		[
			"NAME:Geometry3DCmdTab",
			[
				"NAME:PropServers", 
				"Signal0_2:ThickenSheet:1", 
				"Signal0_1:ThickenSheet:1", 
				"Signal4_3:ThickenSheet:1"
			],
			[
				"NAME:ChangedProps",
				[
					"NAME:Thickness",
					"Value:="		, "-25mm"
				],
				[
					"NAME:Thickness",
					"Value:="		, "-22nm"
				]
			]
		]
	])
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.CreateBox(
	[
		"NAME:BoxParameters",
		"XPosition:="		, xpos,
		"YPosition:="		, ypos,
		"ZPosition:="		, "0mm",
		"XSize:="		, xsize,
		"YSize:="		, ysize,
		"ZSize:="		, "-0.1mm"
	], 
	[
		"NAME:Attributes",
		"Name:="		, "Box2",
		"Flags:="		, "",
		"Color:="		, "(143 175 143)",
		"Transparency:="	, 0,
		"PartCoordinateSystem:=", "Global",
		"UDMId:="		, "",
		"MaterialValue:="	, "\"copper\"",
		"SurfaceMaterialValue:=", "\"\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "0mm",
		"ReferenceTemperature:=", "20cel",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.ChangeProperty(
	[
		"NAME:AllTabs",
		[
			"NAME:Geometry3DCmdTab",
			[
				"NAME:PropServers", 
				"Box2:CreateBox:1"
			],
			[
				"NAME:ChangedProps",
				[
					"NAME:Position",
					"X:="			, xpos,
					"Y:="			, ypos,
					"Z:="			, "-2.2e-05mm"
				]
			]
		]
	])
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.Subtract(
	[
		"NAME:Selections",
		"Blank Parts:="		, "Box1",
		"Tool Parts:="		, "Box2"
	], 
	[
		"NAME:SubtractParameters",
		"KeepOriginals:="	, False,
		"TurnOnNBodyBoolean:="	, True
	])


oEditor.CreateBox(
	[
		"NAME:BoxParameters",
		"XPosition:="		, xpos,
		"YPosition:="		, ypos,
		"ZPosition:="		, "-22nm",
		"XSize:="		, xsize,
		"YSize:="		, ysize,
		"ZSize:="		, "-16nm"
	], 
	[
		"NAME:Attributes",
		"Name:="		, "Box9",
		"Flags:="		, "",
		"Color:="		, "(143 175 143)",
		"Transparency:="	, 0,
		"PartCoordinateSystem:=", "Global",
		"UDMId:="		, "",
		"MaterialValue:="	, "\"copper\"",
		"SurfaceMaterialValue:=", "\"\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "0mm",
		"ReferenceTemperature:=", "20cel",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])
oEditor.AssignMaterial(
	[
		"NAME:Selections",
		"AllowRegionDependentPartSelectionForPMLCreation:=", True,
		"AllowRegionSelectionForPMLCreation:=", True,
		"Selections:="		, "Box9"
	], 
	[
		"NAME:Attributes",
		"MaterialValue:="	, "\"Ge\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "nan ",
		"ReferenceTemperature:=", "nan ",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])


oEditor.CreateBox(
	[
		"NAME:BoxParameters",
		"XPosition:="		, xpos,
		"YPosition:="		, ypos,
		"ZPosition:="		, "-3.8e-05mm",
		"XSize:="		, xsize,
		"YSize:="		, ysize,
		"ZSize:="		, "-1.9um"
	], 
	[
		"NAME:Attributes",
		"Name:="		, "Box3",
		"Flags:="		, "",
		"Color:="		, "(143 175 143)",
		"Transparency:="	, 0,
		"PartCoordinateSystem:=", "Global",
		"UDMId:="		, "",
		"MaterialValue:="	, "\"copper\"",
		"SurfaceMaterialValue:=", "\"\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "0mm",
		"ReferenceTemperature:=", "20cel",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])
oEditor.AssignMaterial(
	[
		"NAME:Selections",
		"AllowRegionDependentPartSelectionForPMLCreation:=", True,
		"AllowRegionSelectionForPMLCreation:=", True,
		"Selections:="		, "Box3"
	], 
	[
		"NAME:Attributes",
		"MaterialValue:="	, "\"SiGe\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "nan ",
		"ReferenceTemperature:=", "nan ",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])

oEditor.CreateBox(
	[
		"NAME:BoxParameters",
		"XPosition:="		, xpos,
		"YPosition:="		, ypos,
		"ZPosition:="		, "-0.001938mm",
		"XSize:="		, xsize,
		"YSize:="		, ysize,
		"ZSize:="		, "-1.4um"
	], 
	[
		"NAME:Attributes",
		"Name:="		, "Box4",
		"Flags:="		, "",
		"Color:="		, "(143 175 143)",
		"Transparency:="	, 0,
		"PartCoordinateSystem:=", "Global",
		"UDMId:="		, "",
		"MaterialValue:="	, "\"copper\"",
		"SurfaceMaterialValue:=", "\"\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "0mm",
		"ReferenceTemperature:=", "20cel",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])
oEditor.AssignMaterial(
	[
		"NAME:Selections",
		"AllowRegionDependentPartSelectionForPMLCreation:=", True,
		"AllowRegionSelectionForPMLCreation:=", True,
		"Selections:="		, "Box4"
	], 
	[
		"NAME:Attributes",
		"MaterialValue:="	, "\"Ge\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "nan ",
		"ReferenceTemperature:=", "nan ",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])

oEditor.CreateBox(
	[
		"NAME:BoxParameters",
		"XPosition:="		, xpos,
		"YPosition:="		, ypos,
		"ZPosition:="		, "-0.003338mm",
		"XSize:="		, xsize,
		"YSize:="		, ysize,
		"ZSize:="		, "-525um"
	], 
	[
		"NAME:Attributes",
		"Name:="		, "Box5",
		"Flags:="		, "",
		"Color:="		, "(143 175 143)",
		"Transparency:="	, 0,
		"PartCoordinateSystem:=", "Global",
		"UDMId:="		, "",
		"MaterialValue:="	, "\"copper\"",
		"SurfaceMaterialValue:=", "\"\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "0mm",
		"ReferenceTemperature:=", "20cel",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])
oEditor.AssignMaterial(
	[
		"NAME:Selections",
		"AllowRegionDependentPartSelectionForPMLCreation:=", True,
		"AllowRegionSelectionForPMLCreation:=", True,
		"Selections:="		, "Box5"
	], 
	[
		"NAME:Attributes",
		"MaterialValue:="	, "\"silicon\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "nan ",
		"ReferenceTemperature:=", "nan ",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])

oEditor.CreateBox(
	[
		"NAME:BoxParameters",
		"XPosition:="		, xpos,
		"YPosition:="		, ypos,
		"ZPosition:="		, "0mm",
		"XSize:="		, xsize,
        "YSize:="		, ysize,
		"ZSize:="		, "900um"
	], 
	[
		"NAME:Attributes",
		"Name:="		, "Box6",
		"Flags:="		, "",
		"Color:="		, "(143 175 143)",
		"Transparency:="	, 0,
		"PartCoordinateSystem:=", "Global",
		"UDMId:="		, "",
		"MaterialValue:="	, "\"copper\"",
		"SurfaceMaterialValue:=", "\"\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "0mm",
		"ReferenceTemperature:=", "20cel",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])
oEditor.AssignMaterial(
	[
		"NAME:Selections",
		"AllowRegionDependentPartSelectionForPMLCreation:=", True,
		"AllowRegionSelectionForPMLCreation:=", True,
		"Selections:="		, "Box6"
	], 
	[
		"NAME:Attributes",
		"MaterialValue:="	, "\"vacuum\"",
		"SolveInside:="		, False,
		"ShellElement:="	, False,
		"ShellElementThickness:=", "nan ",
		"ReferenceTemperature:=", "nan ",
		"IsMaterialEditable:="	, True,
		"IsSurfaceMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,
		"IsLightweight:="	, False
	])