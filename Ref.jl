using LaMEM, GeophysicalModelGenerator, Plots

## define varibles
# slab breakoff depth
breakoff_depth = -100

# Arc crust: 
# Moho temperature
TarcMoho = 700               #  hot arc Moho temperature at 700 C for  hot case, linear gradient of 20C/km


# define crustal thickness
Indian_crust_thickess = 25
Indian_crust_thickess_sediment = 2
Indian_crust_thickess_upper = 12

Indian_margin_crust_thickess = 25
Indian_margin_crust_thickess_sediment = 2
Indian_margin_crust_thickess_upper = 12

Asian_crust_thickess = 35
Asian_crust_thickess_upper = 18

arc_crust_thickness = 40 # km
arc_crust_thickness_upper = 20 #km

reduction_ratio = 0.95 # for plastic softening

# set up model size, horizontal 4000 km and depth to 1200 km 
# also include 10 km air

########################
#  change air thickness from 20 to 10 km
########################
#  May increase local resolution to 1 km (not done yet)
########################
#  need to add passive tracers (not done yet)
########################

model = Model(   Grid(coord_x=[-2000, -200, 200, 2000], bias_x=[0.125,1.0,8.0], nel_x=[200,400,200],
                      coord_y=[-1,1],                   bias_y=[1.0],         nel_y=[1],
                      coord_z=[-1200, -220, 20],        bias_z=[0.2,1.0],     nel_z=[80,240]),

                    PassiveTracers(
                      Passive_Tracer=1, # acticve passive tracer
                      PassiveTracer_Box=[0,600,-1,1,-150,0],  # range of passive tracer coordinate: x_left, x_right, y_negative, y_positive, z_bottom, Z_top
                      PassiveTracer_Resolution  =  [600,1,150]), # resolution in x, y, z direction, the limit of passive tracer is 1e6
    
                    BoundaryConditions( temp_bot        = 1744.0,  # we use 1300C for mantle temp at LAB (at - 90km),  with 0.4 C/km Adiabat gradient, the bottom temp is 1744 C
                                        temp_top        = 20.0,
                                        open_top_bound  = 1),
                    Scaling(GEO_units(  temperature     = 1000,
                                        stress          = 1e9,
                                        length          = 1,
                                        viscosity       = 1e20) )
                     )  

model.Grid

###################
# Non-default setup suggested by Nico
###################

model.ModelSetup = ModelSetup(
                                advect         = "rk2",              # advection scheme
                                interp         = "stag",             # velocity interpolation scheme
                                mark_ctrl      = "subgrid",          # marker control type
                                nmark_lim      = [27, 64],           # min/max number per cell
                                nmark_sub      = 3                   # max number of same phase markers per subcell (subgrid marker control)
                            )

###################
# Model runtime to 20 Myr
# Could be shorter for shallow breakoff
###################

model.Time = Time(  time_end  = 20.0,
                    dt        = 0.0005,  # used to be 0.001 
                    dt_min    = 0.000001,
                    dt_max    = 0.1,
                    nstep_max = 500000,
                    nstep_out = 50,
                    CFL       = 0.2,     # used to be 0.5
                    CFLMAX    = 0.8,
                 )

model.SolutionParams = SolutionParams(  shear_heat_eff = 1.0,  
                                        Adiabatic_Heat = 1.0,
                                        act_temp_diff  = 1,
                                        FSSA = 1.0,     
                                        eta_min   = 5e18,
                                        eta_ref   = 1e21,
                                        eta_max   = 1e25,
                                        min_cohes = 1e3,
                                    )

# To start a test run, we did not include surface processes
# surface level is set to 0 km by default

model.FreeSurface = FreeSurface(    surf_use        = 1,                # free surface activation flag
                                    surf_corr_phase = 1,                # air phase ratio correction flag (due to surface position)
                                    surf_level      = 0.0,              # initial level
                                    surf_air_phase  = 0,                # phase ID of sticky air layer
                                    surf_max_angle  = 40.0 ,             # maximum angle with horizon (smoothed if larger))
                                    sediment_model  = 0,
                                    erosion_model   = 0
                                )

########################
# output strain rate and dev stress tensors
########################
# need to output trackers information
########################

model.Output = Output(  out_density         = 1,
                        out_j2_strain_rate  = 1,
                        out_strain_rate     = 1,
                        out_surf            = 1, 	
                        out_surf_pvd        = 1,
                        out_surf_topography = 1,
                        out_j2_dev_stress   = 1,
                        out_dev_stress      = 1,
                        out_pressure        = 1,
                        out_temperature     = 1,
                        out_ptr          =  1,     # activate
                        out_ptr_ID       =  1,     # ID of the passive tracers
                        out_ptr_phase    =  1,     # phase of the passive tracers
                        out_ptr_Pressure  =  1,     # interpolated pressure
                        out_ptr_Temperature  =  1,     # temperature 
                          )

# define some constant temperatures
# define adiabat mantle temperature gradient : 0.4 C/km

Tair            = 20.0;
Tmantle         = 1300.0;
Adiabat         = 0.4

# set all model temp to mantle temp, will be overwritten later

model.Grid.Temp     .= Tmantle ;            # set mantle temperature (without adiabat at first)
model.Grid.Phases   .= 1;                   # Set Phases to 0 everywhere (#1 is  asthenosphere in this setup)

# phase above Z> 0 km is assigned for air temperature
# phase above Z> 0 km is assigned for air phase (#0 phase)

model.Grid.Temp[model.Grid.Grid.Z .> 0]    .=  Tair;
model.Grid.Phases[model.Grid.Grid.Z .> 0.0 ] .= 0;

# Indian continental plate (major part)
# Crustal thickness (35) and LAB depth (135) from Wang et al., 2023. Secular craton evolution due to cyclic deformation of underlying dense mantle lithosphere. Nature Geoscience.

add_box!(model;  xlim    = (-2000.0, -750.0), 
                ylim    = (model.Grid.coord_y...,), 
                zlim    = (-120.0, 0.0),
                Origin  = nothing, StrikeAngle=0, DipAngle=0,
                phase   = LithosphericPhases(Layers=[Indian_crust_thickess_sediment Indian_crust_thickess_upper Indian_crust_thickess-Indian_crust_thickess_sediment-Indian_crust_thickess_upper 120-Indian_crust_thickess], Phases=[14 15 16 3 1] ),
                T       = LinearTemp(               Ttop        = Tair,
                                                    Tbot        = Tmantle
                                                     ))

# Indian continental plate (thinned part)


add_box!(model;  xlim    = (-750.0, -250.0), 
ylim    = (model.Grid.coord_y...,), 
zlim    = (-120, 0.0),
Origin  = nothing, StrikeAngle=0, DipAngle=0,
phase   = LithosphericPhases(Layers=[Indian_margin_crust_thickess_sediment Indian_margin_crust_thickess_upper Indian_margin_crust_thickess-Indian_margin_crust_thickess_sediment-Indian_margin_crust_thickess_upper 120-Indian_margin_crust_thickess], Phases=[14 15 16 3 1] ),
T       = LinearTemp(                Ttop        = Tair,
                                     Tbot        = Tmantle   
                                     ) )                                                  

# Asian continental plate (thinned part)
# Slightly thinned compared to the rest of the Asian plate (30 crust + 70 mantle lithosphere)

add_box!(model;  xlim    = (250.0, 350.0), 
                ylim    = (model.Grid.coord_y...,), 
                zlim    = (-100, 0.0),
                Origin  = nothing, StrikeAngle=0, DipAngle=0,
                phase   = LithosphericPhases(Layers=[15 15 70], Phases=[12 13 4 1] ),   
                T       = LinearTemp(                Ttop        = Tair,
                                                     Tbot        = Tmantle  
                                                     ) )

# Asian continental plate (the major part)
# Crustal thickness (35) and mantle lithosphere (100) from Wang et al., 2023. Secular craton evolution due to cyclic deformation of underlying dense mantle lithosphere. Nature Geoscience.

add_box!(model;  xlim    = (350.0, 2000.0), 
                ylim    = (model.Grid.coord_y...,), 
                zlim    = (-135, 0.0),
                Origin  = nothing, StrikeAngle=0, DipAngle=0,
                phase   = LithosphericPhases(Layers=[Asian_crust_thickess_upper Asian_crust_thickess-Asian_crust_thickess_upper 135-Asian_crust_thickess], Phases=[12 13 4 1] ),
                T       = LinearTemp(                Ttop        = Tair,
                                                     Tbot        = Tmantle
                                                     ) )


# Arc phase, from trench 200-300 km
# trench is 100 km left of the leftmost margin of Asian plate due to accretionary complex
# so arc is 100-200 km from the leftmost margin (x=250) of Asian plate (so arc location is x=350-450)

###############################################################################
# Crust and mantle lithosphere are added seperately to set different temperatures
###############################################################################

# Arc crust: 


add_box!(model;  xlim    = (350.0, 450.0), 
                ylim    = (model.Grid.coord_y...,),  
                zlim    = (-arc_crust_thickness, 0.0),                          # from 0 to -35 km, 35 km thick
                Origin  = nothing, StrikeAngle=0, DipAngle=0,
                phase   = LithosphericPhases(Layers=[arc_crust_thickness_upper arc_crust_thickness-arc_crust_thickness_upper], Phases=[17 18 8] ),                  # arc phase #15, 16, 17
                T       = LinearTemp(Ttop        = Tair,
                                     Tbot        = TarcMoho
                ))

# Arc mantle lithophere: 

add_box!(model;  xlim    = (350.0, 450.0), 
                ylim    = (model.Grid.coord_y...,), 
                zlim    = (-arc_crust_thickness- 100, -arc_crust_thickness),                         # from - 35 to -135 km, 100 km thick
                Origin  = nothing, StrikeAngle=0, DipAngle=0,
                phase   = ConstantPhase(8),                    # arc mantle lithophere phase #17
                T       = LinearTemp(Ttop        = TarcMoho,
                                     Tbot        = Tmantle
                ))

# This is oceanic plate to be convergened (500 km wide)
# 90 km thick for the whole lithosphere consistant with halfspace cooling age at ~80 Ma

add_box!(model;  xlim    = (-250.0, 250.0), 
                ylim    = (model.Grid.coord_y..., ), 
                zlim    = (-90.0, 0.0),
                Origin  = nothing, StrikeAngle=0, DipAngle=0,
                phase   = LithosphericPhases(Layers=[7 83], Phases=[10 2 1] ),
                T       = HalfspaceCoolingTemp(    Tsurface    = Tair,
                                                   Tmantle     = Tmantle,
                                                   Age         = 80 ))

# Weak zone on top of pre-subducted oceanic plate to faciliate initial decoupling

add_box!(model;  xlim    = (260, 550), 
                ylim    = (model.Grid.coord_y...,), 
                zlim    = (-5.0, 0.0),
                Origin  = nothing, StrikeAngle=0, DipAngle=30,
                phase   = ConstantPhase(6),
                T       = nothing ) 

# the starting position is at x= 250 km, given the 30 deg angle and the top is at 130 km depth

add_box!(model;  xlim    = (250, 550), 
                ylim    = (model.Grid.coord_y...,), 
                zlim    = (-90.0, 0.0),
                Origin  = nothing, StrikeAngle=0, DipAngle=30,
                phase   = LithosphericPhases(Layers=[7 83], Phases=[10 2 1], Tlab=Tmantle ),
                T       = HalfspaceCoolingTemp(     Tsurface    = Tair,
                                                    Tmantle     = Tmantle,
                                                    Age         = 80      ) )

# control box: oceanic crust part # 12 top 7 km 

add_box!(model;  xlim    = (-250, -240), 
                ylim    = (model.Grid.coord_y...,), 
                zlim    = (-7.0, 0.0),
                Origin  = nothing, StrikeAngle=0, DipAngle=0,
                phase   = ConstantPhase(11),
                T       = nothing ) 

# control box: oceanic lithospheric mantle part # 13

add_box!(model;  xlim    = (-250, -240), 
                ylim    = (model.Grid.coord_y...,), 
                zlim    = (-90.0, -7.0),
                Origin  = nothing, StrikeAngle=0, DipAngle=0,
                phase   = ConstantPhase(5),
                T       = nothing ) 


# add Adiabat temperture profile when Z < -90 km, use the LAB depth of oceanic lithosphere

model.Grid.Temp[model.Grid.Grid.Z .< -90]   .=  model.Grid.Temp[model.Grid.Grid.Z.< -90] .- model.Grid.Grid.Z[model.Grid.Grid.Z.< -90].*Adiabat;  # "-" sign is used since Z is a negative value(corrdinate)

plot_cross_section(model, y=0, field=:phase)


plot_cross_section(model, y=0, field=:temperature)


# Softening
softening =       Softening(    ID   = 0,   			# softening law ID
                                APS1 = 0.1, 			# begin of softening APS
                                APS2 = 0.5, 			# end of softening APS
                                A    = reduction_ratio, 		    # reduction ratio
)

air         = Phase(    Name        = "air",                                     
                        ID          = 0,                                                # phase id 
                        rho         = 50.0,                                             # density [kg/m3]                                           # coeff. of thermal expansion [1/K]
                        eta         = 1e19,
                        G           = 5e10,                                             # elastic shear module [Pa]
                        k           = 100,                                              # conductivity
                        Cp          = 1e6,                                              # heat capacity
                        ch          = 10e6,                                             # cohesion [MPa]
                        fr          = 0.0,                                              # friction angle	
                    )

dryPeridotite = Phase(  Name        = "dryPeridotite",                                     
                        ID          = 1,                                                # phase id  [-]
                        rho         = 3300.0,                                           # density [kg/m3]
                        alpha       = 3e-5,                                             # coeff. of thermal expansion [1/K]
                        disl_prof   = "Dry_Olivine_disl_creep-Hirth_Kohlstedt_2003",
                        #Vn          = 14.5e-6,
                        diff_prof   = "Dry_Olivine_diff_creep-Hirth_Kohlstedt_2003",
                        #Vd          = 14.5e-6,  
                        #peir_prof   = "Olivine_Peierls-Kameyama_1999",                  # PEIERLS creep profile
                        G           = 5e10,                                             # elastic shear module [Pa]
                        k           = 3,                                                # conductivity
                        Cp          = 1000.0,                                           # heat capacity
                        ch          = 30e6,                                             # cohesion [Pa]
                        fr          = 20.0,                                             # friction angle	
                        A           = 6.6667e-12,                                       # radiogenic heat production [W/kg], which is 0.02 uW/m3
                        chSoftID    = 0,      	                                        # cohesion softening law ID
                        frSoftID    = 0,      	                                        # friction softening law ID
                        )

oceanicCrust = Phase(   Name        = "oceanCrust",                                     
                        ID          = 10,                                                # phase id  [-]
                        rho         = 3000.0,                                           # density [kg/m3]
                        alpha       = 3e-5,                                             # coeff. of thermal expansion [1/K]
                        disl_prof   = "Plagioclase_An75-Ranalli_1995",
                        G           = 5e10,                                             # elastic shear module [Pa]
                        k           = 3,                                                # conductivity
                        Cp          = 1000.0,                                           # heat capacity
                        ch          = 5e6,                                              # cohesion [Pa]
                        fr          = 0.0,                                              # friction angle	
                        A           = 6.6667e-11,                                       # radiogenic heat production [W/kg], which is 0.2 uW/m3
                     )

oceanicLithosphere = copy_phase(    dryPeridotite,
                                    Name            = "oceanicLithosphere",
                                    ID              = 2
                                )

AsianContinentalCrust_upper = copy_phase(      oceanicCrust,
                                    Name            = "AsianContinentalCrust_upper",
                                    ID              = 12,
                                    disl_prof       = "Quarzite-Ranalli_1995",
                                    rho             = 2800.0,  
                                    ch              = 30e6,
                                    fr              = 20.0,
                                    A         	    = 5.3571e-10,  #radiogenic heat production [W/kg], which is 1.5 uW/m3
                                    chSoftID  	    = 0,      	                                        
                                    frSoftID  	    = 0,      	
                                    
                             )

AsianContinentalCrust_lower = copy_phase(      oceanicCrust,
                                    Name            = "AsianContinentalCrust_lower",
                                    ID              = 13,
                                    disl_prof       = "Quarzite-Ranalli_1995",
                                    rho             = 2800.0,  
                                    ch              = 30e6,
                                    fr              = 20.0,
                                    A         	    = 5.3571e-10,  #radiogenic heat production [W/kg], which is 1.5 uW/m3
                                    chSoftID  	    = 0,      	                                        
                                    frSoftID  	    = 0,  
                                    
                             )

AsianContinentalMantleLithosphere = copy_phase(    dryPeridotite,
                                    Name            = "AsianContinentalMantleLithosphere",
                                    ID              = 4
                                )

IndianContinentalCrust_sediment = copy_phase(      oceanicCrust,
                                    Name            = "IndianContinentalCrust_sediment",
                                    ID              = 14,
                                    disl_prof       = "Quarzite-Ranalli_1995",
                                    rho             = 2700.0,  
                                    ch              = round(30e6 * (1 - reduction_ratio); digits=1),
                                    fr              = round(20.0 * (1 - reduction_ratio); digits=1),
                                    A         	    = 5.3571e-10,  #radiogenic heat production [W/kg], which is 1.5 uW/m3
                                   #  chSoftID  	    = 0,      	                                        
                                   #  frSoftID  	    = 0,   
                                    
                             )

IndianContinentalCrust_upper = copy_phase(      oceanicCrust,
                                    Name            = "IndianContinentalCrust_upper",
                                    ID              = 15,
                                    disl_prof       = "Quarzite-Ranalli_1995",
                                    rho             = 2700.0,  
                                    ch              = 30e6,
                                    fr              = 20.0,
                                    A         	    = 5.3571e-10,  #radiogenic heat production [W/kg], which is 1.5 uW/m3
                                    chSoftID  	    = 0,      	                                        
                                    frSoftID  	    = 0,   
                                    
                             )

IndianContinentalCrust_lower = copy_phase(      oceanicCrust,
                                    Name            = "IndianContinentalCrust_lower",
                                    ID              = 16,
                                    disl_prof       = "Plagioclase_An75-Ranalli_1995",
                                    rho             = 2900.0,  
                                    ch              = 30e6,
                                    fr              = 20.0,
                                    A         	    = 5.3571e-10,  #radiogenic heat production [W/kg], which is 1.5 uW/m3
                                    chSoftID  	    = 0,      	                                        
                                    frSoftID  	    = 0,      
                                                                            
                             )

IndianContinentalMantleLithosphere = copy_phase(    dryPeridotite,
                                    Name            = "IndianContinentalMantleLithosphere",
                                    ID              = 3
                                )

eclogite      = Phase(  Name        = "eclogite",                                     
                        ID          = 9,                                                  # phase id  [-]

                        # some testing rho and eta 
                        rho         = 3500.0,                                             # density [kg/m3]                                           # coeff. of thermal expansion [1/K]
                        eta         = 1e22,

                        # same properties as oceanic crust
                        G           = 5e10,                                             # elastic shear module [Pa]
                        k           = 3,                                                # conductivity
                        Cp          = 1000.0,                                           # heat capacity
                        ch          = 5e6,                                              # cohesion [Pa]
                        fr          = 0.0,                                              # friction angle	
                        A           = 6.6667e-11,                                        # radiogenic heat production [W/kg]
                    )

weakzone_slab     = Phase(  Name        = "weakzone_slab",                                     
                        ID          = 6,                                                # phase id  [-]

                        # some testing rho and eta 
                        rho         = 3000,                                             # density [kg/m3]                                           # coeff. of thermal expansion [1/K]
                        eta         = 1e19,

                        # same properties as oceanic crust
                        G           = 5e10,                                             # elastic shear module [Pa]
                        k           = 3,                                                # conductivity
                        Cp          = 1000.0,                                           # heat capacity
                        ch          = 10e6,                                             # cohesion [Pa]
                        fr          = 0.0,                                              # friction angle	
                        A           = 0,                                                # radiogenic heat production [W/kg]
                    )

# control box same properties as oceanic crust. oceanic_crust
controlbox_crust   =  copy_phase(    oceanicCrust,
                                     Name            = "controlbox_crust",
                                     ID              = 11
                                )

# control box same properties as oceanic mantle lithosphere 
controlbox_mantle   =  copy_phase(    dryPeridotite,
                                      Name            = "controlbox_mantle",
                                      ID              = 5
                                  )

weakzone_box     = Phase(  Name        = "weakzone_box",          # same as mantle                             
                           ID          = 7,                                               # phase id  [-]
                            # some testing rho and eta 
                           rho         = 3300,                                             # density [kg/m3]                                           # coeff. of thermal expansion [1/K]
                           eta         = 1e19,
                           alpha       = 3e-5,
                           # same properties as manlte phase
                           G           = 5e10,                                             # elastic shear module [Pa]
                           k           = 3,                                                # conductivity
                           Cp          = 1000.0,                                           # heat capacity
                           ch          = 30e6,                                             # cohesion [Pa]
                           fr          = 20.0,                                             # friction angle	
                           A           = 6.6667e-12,                                       # radiogenic heat production [W/kg]
                          )

# here is added arc crust
arcCrust_upper = copy_phase(      oceanicCrust,
                                    Name            = "arcCrust_upper",
                                    ID              = 17,
                                    disl_prof       = "Quarzite-Ranalli_1995",
                                    rho             = 2800.0,  
                                    ch              = 30e6,
                                    fr              = 20.0, 
                                    A         	    = 5.3571e-10,
                                    chSoftID  	    = 0,      	                                        
                                    frSoftID  	    = 0,   
                                                                          
                             )

# here is added arc crust
arcCrust_lower = copy_phase(      oceanicCrust,
                                    Name            = "arcCrust_lower",
                                    ID              = 18,
                                    disl_prof       = "Quarzite-Ranalli_1995",
                                    rho             = 2800.0,  
                                    ch              = 30e6,
                                    fr              = 20.0, 
                                    A         	    = 5.3571e-10,
                                    chSoftID  	    = 0,      	                                        
                                    frSoftID  	    = 0,      	   
                                                                  
                             )

ArcMantleLithophere = copy_phase(    dryPeridotite,
                                        Name            = "ArcMantleLithophere",
                                        ID              = 8
                                )

# Add phase transitions
# Depth dependent phase transition, oceanic crust #9  -> eclogite #10
PT0 = PhaseTransition(ID=0, Type="Constant", Parameter_transition="Depth",  PhaseBelow = [9], PhaseAbove=[10], PhaseDirection="AboveToBelow", ConstantValue=-80, ResetParam="APS")

# Depth dependent phase transition, breakoff control box, oceanic crust -> weak zone_box 
PT1 = PhaseTransition(ID=1, Type="Constant", Parameter_transition="Depth",  PhaseBelow = [7], PhaseAbove=[11], PhaseDirection="AboveToBelow", ConstantValue= breakoff_depth, ResetParam="APS")

# Depth dependent phase transition, breakoff control box, oceanic lithosphere mantle -> weak zone_box
PT2 = PhaseTransition(ID=2, Type="Constant", Parameter_transition="Depth",  PhaseBelow = [7], PhaseAbove=[5], PhaseDirection="AboveToBelow", ConstantValue= breakoff_depth, ResetParam="APS")

model.Materials.PhaseTransitions = [PT0, PT1, PT2]

rm_phase!(model)
add_phase!( model, 
            air,
            dryPeridotite,
            oceanicCrust,
            oceanicLithosphere,
            AsianContinentalCrust_upper,
            AsianContinentalCrust_lower,
            AsianContinentalMantleLithosphere,
            IndianContinentalCrust_sediment,
            IndianContinentalCrust_upper,
            IndianContinentalCrust_lower,
            IndianContinentalMantleLithosphere,
            eclogite,
            weakzone_slab,
            weakzone_box,
            controlbox_crust,
            controlbox_mantle,
            arcCrust_upper,
            arcCrust_lower,
            ArcMantleLithophere
          )

add_softening!( model, softening)

model.Solver = Solver(  SolverType      = "multigrid",
                        MGLevels        = 3,
                        MGCoarseSolver 	= "mumps",
                        DirectSolver 	= "mumps",
                        PETSc_options   = [ "-snes_ksp_ew",
                                            "-snes_ksp_ew_rtolmax 1e-4",
                                            "-snes_rtol 5e-3",			
                                            "-snes_atol 1e-4",
                                            "-snes_max_it 200",
                                            "-snes_PicardSwitchToNewton_rtol 1e-3", 
                                            "-snes_NewtonSwitchToPicard_it 20",
                                            "-js_ksp_type fgmres",
                                            "-js_ksp_max_it 20",
                                            "-js_ksp_atol 1e-8",
                                            "-js_ksp_rtol 1e-4",
                                            "-snes_linesearch_type l2",
                                            "-snes_linesearch_maxstep 10",
                                            "-da_refine_y 1"
                                        ]
                    )

prepare_lamem(model, 8);

# ── Append pushing-block boundary conditions to output.dat ───────────────────

sleep(1)  # ensure output.dat is fully written before reading

original_file = "output.dat"

pushing_block_text = """

# ----------------------------------------------------
# Bezier/Pushing block
# ----------------------------------------------------

# velocity    - cm/yr , 5cm/yr = 50 mm/yr  = 50 km/Myr

#  2 stage velocity, push from both sides:
#
#  0-5 Myr: 10 cm/yr total, 5 cm/yr (50 km/Myr) each side
#  initial touch of two continents at 5 Myr
#  5-10 Myr: 5 cm/yr total, 2.5 cm/yr (25 km/Myr) each side


<BCBlockStart>

# Right side push block
npath		        	=   3
time                    =   0 5 20
theta                   =   0 0 0
path                    =   2000 0 1750 0 1375 0
# polygon
npoly                   =   4
poly                    =   1900 -1 2000 -1 2000 1 1900 1
bot                     =   -135.0
top                     =   0.0
<BCBlockEnd>

# Left side push block

<BCBlockStart>
npath			        =   3
time                    =   0 5 20
theta                   =   0 0 0
path                    =   -2000 0 -1750 0 -1375 0
# polygon
npoly                   =   4
poly                    =   -2000 -1 -1900 -1 -1900 1 -2000 1
bot                     =   -135.0
top                     =   0.0
<BCBlockEnd>

"""

original_content = read(original_file, String)
write(original_file, original_content * pushing_block_text)

println("Pushing block appended to '$original_file'")

run_lamem(model, 8)