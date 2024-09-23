
using Base.Threads
using UnicodePlots
using Distributions
using LinearAlgebra
using StatsPlots
using JSON
using JSON3
#include("functions_list_herb.jl")
#using .functions_list_herb
########################################## FUNÇÕES ALOMÉTRICAS #############################

println(nthreads())


function find_velocity(mass)
    # mass in kg
    # from kramer 2010 "allometric scaling of resource acquisition"
    #Consumer Velocity (meters/second)
    velocity = (0.5 * mass^0.13);
    return velocity
end

function bite_size_allo(mass)
    # mass in kg
    # bite size in kg
    bite_size = (0.002 * (mass)^0.969)/1000; # [kg]
    return bite_size

end


function bite_rate_allo(mass)
    # Why elephants have trunks Pretorius 2015
    bite_rate = 0.37 * mass^(-0.024) #(bites/s)  
    return bite_rate
end


function alpha_allo(mass)
    bite_rate = bite_rate_allo(mass); # mass in kg, bite/s
    bite_size = bite_size_allo(mass) # mass in kg, bite size in kg
    alpha = bite_rate * bite_size # bite/s * kg/bite
    return alpha
end



function number_of_chews(mass)
    # shipley 94
   #chew/g (processed to mean particle size (allo) and swallowed)
   # mass in kg
    chews_per_gram = 343.71* 1000 * (mass^(-0.83)); 
    return chews_per_gram
 end

 function chew_rate_allo(mass)
    # from "dental functional morphology predicts scaling"
    # mass in kg
    # duration in ms
     chewing_cycle_duration = (228.0* (mass)^0.246) / 1000;
     return 1 / (chewing_cycle_duration )  #[s/chew -> chews/s]
 
 end


 function chew_allo(mass)
    # allometric function for mouth->gut rate
   chew_rate = chew_rate_allo(mass); # chews/s
   chews_per_gram = number_of_chews(mass)   # chew/kg
   chew = chew_rate / chews_per_gram        # chew/s / chew/kg -> kg/s
   return chew
end


function mean_retention_time(mass)
    # from Muller 2013 - Assessing the Jarmain-Bell Principle
    # mrt in hrs
    # bm in kg
    # mean retention of a particle in the gut [s]

    mean_retention_time = (30.3 * (mass)^0.109)
    return mean_retention_time * 60 * 60 # [hr -> s]

end

function gut_volume_g(mass)
    capacity = (0.030 * (mass)^0.881);
    #return capacity * 1000 #[kg -> g]
    return capacity
end



function mean_particle_mass(mass)
    # from "comparative chewing efficiency"
    # mean particle size in mm
    # mass in g

    mass_g = mass * 1000; # convert mass kg -> g
    mean_particle_size = 0.0;
    mean_particle_size = (6.61 * (mass_g)^0.26)
    #mean_particle_size = (6.61 * (mass)^0.26)
    volume = (4/3) * pi * (1/2 * mean_particle_size)^3; # [mm^3]
    particle_mass = 0.0004 * volume; # [g/mm^3 * mm^3 = g]

    return particle_mass

end
#####################################################################################################

nr=2

max_m = 1; # Maximum mean abundance
min_m = 0.1; # Minimum mean abundance

m = collect(min_m:(max_m-min_m)/(nr-1):max_m)

v = 0.00001

alpha = m.^2 ./v

mdist = alpha ./ (m .* (alpha .- 1)); # Mean distance to resources

max_res_mcal= 1;
min_res_mcal= 0.5;

res_mcal = reverse(collect(max_res_mcal:-(max_res_mcal-min_res_mcal)/(length(m)-1):min_res_mcal)) #quantidade de kj de cada recurso
targetvalues = collect(0.75:0.25:1.0); # Possible targeting values
tweight = [0; repeat(targetvalues, outer=(nr,1))]; # Targeting weights
tid = [0; repeat(collect(1:nr), inner=(length(targetvalues),1))]; # Resource IDs for targeting
tinfo = Tuple([tid, tweight]); # não serve para nada
#println(tinfo)
activehours = 5; # Foraging time in hours
tmax_bout = activehours * 60 * 60; # Convert foraging time to seconds

#################################### DAILY SIMULATION ######################################
function dailysim(nr,alpha,m,res_mcal,mass,tmax_bout,tid,tweight,target,debug)

    velocity = find_velocity(mass) # Forager's velocity in m/s
    beta = bite_size_allo(mass);
    edensity = 0.0
    # mouthrate = bite_rate * bite_size; # bite/s * g/bite = grams/s
    # 1/mouthrate is seconds/1 gram

    chewrate = chew_allo(mass); #g/s
    tchewgram = 1/chewrate; #s/g
    tchew = tchewgram * beta; #s/g * g/bite = s/bite

    
    alpha = Array(alpha);
    # c = Array(c);
    m = Array(m);
    #Negative binomial parameters
    # r = alpha;
    # p = c ./ (c .+ 1);

    gammadist = Gamma.(alpha, m ./ alpha);

    t_travel = Array{Float64}(undef,1);
    t_chew = Array{Float64}(undef,1);
    bites = Array{Float64}(undef,1);
    GUT = Array{Float64}(undef,1);
    # propres = Array{Float64}(undef,nr);
    max_gut = gut_volume_g(mass)
    bernoulidist = Bernoulli(1-tweight[target]); #prob of targeting the closest food or the chosed food
            
    #data = zeros(Float64,configurations);
    
    #Slows down the consumer that is generalizing on foods (search)
    modvelocity = maximum([tweight[target],1/nr])*velocity; #slowering down the generalists
    
    t_travel[1] = 0.0;
    t_chew[1] = 0.0;
    bites[1] = 0.0;
    
    number_of_successes = zeros(Int64,nr); #number of encounters with each iten
    nearest_resource = 0;
    t=0.0;
    distance_to_resource = zeros(Float64,nr);
    nearest_distance = 0.0;    
    gut = 0.0
    #WATTS/KG
    b0_bmr = 0.018; # Basal metabolic rate constant (watts g^-0.75)
    b0_fmr = 0.047; # Field metabolic rate constant (watts g^-0.75)
    basal_mr= (b0_bmr*(1000^0.75))*mass^0.75; # Basal metabolic rate for kg
    #basal_mr = (b0_bmr * (mass^0.75)); # Basal metabolic rate
    field_mr = (b0_fmr*(1000^0.75))*mass^0.75; # Active metabolic rate for kg
    #cost_whr = (b0_fmr*(mass^0.75))*activehours + (b0_bmr*(mass^0.75))*nonactivehours; #watt*hour
    #cost_mcal=cost_whr*8.598*10^-4
    
    #WATTS/KG
    # b0_bmr = 18#/4184000; # Basal metabolic rate constant (Mcal/s Kg^-0.75)
    # b0_fmr = 47#/4184000; # Field metabolic rate constant (Mcal/s Kg^-0.75)
    # basal_mr = (b0_bmr * (mass^0.75)); # Basal metabolic rate
    # field_mr = (b0_fmr * (mass^0.75)); # Active metabolic rate
    #Metabolic cost in watt*seconds
    cost_ws = 0.0;


    while t < tmax_bout
    
        for i=1:nr
            
            distance_to_resource[i] = rand(Exponential(1.0/rand(gammadist[i])));
            
        end
    
        # if debug==true
        #     println("vetor de distancia até recurso ", distance_to_resource)
        # end
        #calculate the distance for every resource
        #ALTERNATIVE
        distancetuple = findmin(distance_to_resource); #finding the closest resource and resource type
        nearest_distance = distancetuple[1];
        nearest_resource = distancetuple[2];
        # if debug==true 
        #     println("dist min: ",distancetuple[1])
        #     println("comida mais perto: ", distancetuple[2])
        # end
        
        if rand(bernoulidist) == 0
    
            # if debug==true
            #     println("Huumm... bolo de murango")
            # end
            #The rodent will move towards the Closest of targeted resources regardless if it's the closest as a whole
            specializedtarget = tid[target]; #********* target tem que ter tamanho 101
        
            sdistance_to_resource = distance_to_resource[specializedtarget]; #******** mas distance_to_resource tem tamanho 10
    
            # if debug==true
            #     println("sdistance_to_resource ", sdistance_to_resource)
            # end
            
            sdistance_traveled = sum(sdistance_to_resource);
    
            deltat = sdistance_traveled/modvelocity;
            t += deltat;
            t_travel[1] += deltat;
    
            
            #Obtains the resource if there is time left in tmax_bout
            if (tmax_bout > (t+ deltat)) && (max_gut > gut)
    
                # if debug == true
                #     println("Oi geeente, vamos aumossar?")
                # end 
                #If not an insect, success is gauranteed
                # if in(7,tid[target]) == false
                number_of_successes[specializedtarget] += 1;
                t += tchew; #time
                t_chew[1] += tchew; #time
    
                # Pass Mouth-Unit to Gut Total (boundary conditions in across-day sim)
                # resgain is kJ per volume (mouth volume)
                gut += beta; #kg/bite
                bites .+ 1;
                cost_ws += field_mr*tchew; #cost of handling time
    
                GUT[1] = gut; 
                edensity += GUT[1]*res_mcal[specializedtarget] #gain energy
    
            else 
                break
                
            end
        
        else
            # if debug==true
            #     println("To cagado de fome")
            # end
            #The rodent will move towards the closest resource
            deltat = nearest_distance/modvelocity;
            t += deltat;
            t_travel[1] += deltat;
            cost_ws += field_mr*deltat;
            
            if (tmax_bout > (t + deltat)) && (max_gut > gut)
                number_of_successes[nearest_resource] += 1;
                t += tchew; #time
                t_chew[1] += tchew; #time
                gut += beta; #grams/bite
                bites .+ 1;
                cost_ws += field_mr*tchew;
                GUT[1] = gut
                edensity += GUT[1]*res_mcal[nearest_resource]   
    
            else
                break
                
            end
    
        end
    end
    total_t = 60*60*24;
    rest_t = total_t - t;
    cost_ws += basal_mr*rest_t;
    #println(rest_t)
    #Convert cost to Mcal
    cost_mcal = cost_ws*8.598*10^-4;

    total_mcal=edensity
    propres = ((res_mcal).*number_of_successes);#mcal per food iten
    return total_mcal,propres,number_of_successes,cost_mcal;
end

########################################################################################################
mass_dict = Dict()
#mass[strat] = Dict()
configurations = 20000
nbins=50
##################################### HISTOGRAM ###################################################

function calculate_histogram(data, nbins)
    # Find the minimum and maximum values in the data to determine the range.
    min_val = minimum(data)
    max_val = maximum(data)
    # Create bin edges that span from the minimum to the maximum value, divided into nbins.
    bin_edges = range(min_val, max_val, length=nbins+1)
    # Initialize a vector to count the number of data points in each bin.
    bin_counts = zeros(Int, nbins)

    # Iterate through each data point and determine which bin it belongs to.
    for val in data
        for i = 1:nbins
            if val >= bin_edges[i] && val < bin_edges[i+1]
                bin_counts[i] += 1
                break
            end
        end
    end

    # Special case: if the last data point equals the maximum value, add it to the last bin.
    if data[end] == max_val
        bin_counts[end] += 1
    end

    return bin_edges, bin_counts
end

############################################# DETERMINE BINS ##################################################

function find_bin_index(value, bins)
    nbins = length(bins) - 1
    # Loop through the bins to find where the value fits.
    for i = 1:nbins
        if value >= bins[i] && value < bins[i+1]
            return i
        end
    end
    # If the value is exactly equal to the maximum, assign it to the last bin.
    return nbins
end

########################################## MAIN LOOP DAILY SIMU ###########################################
function call_dailysimu(m_min, m_max)
    #Threads.@threads for mass in m_min:1:m_max
    Threads.@threads for mass in [m_min, m_max]
        mass_dict[mass] = Dict()
        for target in 1:length(tid) 
            mass_dict[mass][target] = Dict()
            gains = Float64[]
            costs = Float64[]
            for _ in 1:configurations
                # Call the dailysim function to simulate one round of foraging and collect the total kilojoules gained and the cost.
                total_mcal, _, _, cost_mcal = dailysim(nr, alpha, m, res_mcal, mass, tmax_bout, tid, tweight, target, false)
                # Store the results in the gains and costs vectors.
                push!(gains, total_mcal)
                push!(costs, cost_mcal)
            end
            gains_bin_edges, gains_bin_counts = calculate_histogram(gains, nbins)
            costs_bin_edges, costs_bin_counts = calculate_histogram(costs, nbins)
            gains_probabilities = gains_bin_counts / sum(gains_bin_counts)
            costs_probabilities = costs_bin_counts / sum(costs_bin_counts)

            # Calculate the midpoints of bins for gains and costs for better representation.
            gains_bin_midpoints = (gains_bin_edges[1:end-1] + gains_bin_edges[2:end]) / 2;
            costs_bin_midpoints = (costs_bin_edges[1:end-1] + costs_bin_edges[2:end]) / 2;
        
            # Pair up each gain with its corresponding cost.
            gain_cost_pairs = [(gains[i], costs[i]) for i in 1:length(gains)]
            
            # Define bins for the joint distribution of gains and costs.
            gain_bins = range(minimum(gains), maximum(gains), length=nbins+1)
            cost_bins = range(minimum(costs), maximum(costs), length=nbins+1)
            
            # Initialize a matrix to calculate the joint histogram of gains and costs.
            joint_histogram = zeros(Int, nbins, nbins)
            
            # Fill in the joint histogram by finding the appropriate bin for each gain-cost pair.
            for pair in gain_cost_pairs
                gain_bin_index = find_bin_index(pair[1], gain_bins)
                cost_bin_index = find_bin_index(pair[2], cost_bins)
                joint_histogram[gain_bin_index, cost_bin_index] += 1
            end
            
            # Convert the counts in the joint histogram to probabilities.
            joint_probabilities = joint_histogram / sum(joint_histogram)
            mass_dict[mass][target]["gains"]= gain_bins
            mass_dict[mass][target]["costs"]= cost_bins
            mass_dict[mass][target]["prob"]= joint_probabilities
        end
    end
    return mass_dict
end
############################################ SDP MODEL #####################################################

function bc(yvalue, yc, ymax)
    return clamp(yvalue, yc, ymax)
end ####### FUNÇÃO DE CORTE #############

# function interpolate(xvalue, nvalue, W, t)
#     lowx = floor(Int, xvalue)
#     highx = lowx + 1
#     px = highx - xvalue
#     lown = floor(Int, nvalue)
#     highn = lown + 1
#     pn = highn - nvalue

#     #println(lowx)
#     #println(W[lowx,lown,t+1])
#     Winterp = 0.0
#     if pn < 1
#         if px == 1
#             Winterp += px * pn * W[lowx, lown, t+1]
#             Winterp += px * (1-pn) * W[lowx, highn, t+1]
#         else
#             # Winterp += px * pn * W[lowx, lown, t+1]
#             # Winterp += (1-px) * pn * W[highx, lown, t+1]
#             # Winterp += (1-px) * (1-pn) * W[highx, highn, t+1]
#             Winterp += px * pn * W[lowx,lown,t+1]
#             Winterp += px * (1-pn) * W[lowx,highn,t+1]
#             Winterp += (1-px) * pn * W[highx,lown,t+1]
#             Winterp += (1-px) * (1-pn) * W[highx,highn,t+1]
#         end
#     else
#         if px == 1
#             Winterp += px * pn * W[lowx, lown, t+1]
#         else
#             Winterp += px * pn * W[lowx,lown,t+1]
#             Winterp += (1-px) * pn * W[highx,lown,t+1]
#         end
#     end

#     return Winterp
# end
m_min = 700
m_max = 1500
println("início dailysim")
mass_dict = call_dailysimu(m_min, m_max)

println("fim dailysim")
######################################## SDP FUNCTION ##########################################

function SDP(m_min, m_max)
    S = Dict()
    #D = Dict() 
    #R = Dict() 
    T = Dict()# Assuming D holds integer indices
    tmax = 30
    d = 0.001
    xc = 0
    Threads.@threads for m in [m_min, m_max]
    #Threads.@threads for m in m_min:1:m_max
        ymax = Int64(round(4.87*0.02*(m^1.19)))
        yc = Int64(round(4.87*0.1*0.02*(m^1.19)))
        xmax=Int64(round(gut_volume_g(m)*4.87))
        # Maximum gut capacity        # NOTE THIS IS TOO LONG
        S[m] = zeros(Float64, xmax, ymax, tmax)
        #D[m] = zeros(Float64, xmax, ymax, tmax)#probability of a chosen resource
        #R[m] = zeros(Float64, xmax, ymax, tmax)#resource
        T[m] = zeros(Float64, xmax, ymax, tmax)

        #mrt = mean_retention_time(m); # seconds / PARTICLE
        
        #particle_mass = mean_particle_mass(m); #gram / particle
        
        # Passage rate of food (rate of flow from gut to body)
        #passrate = (1/mrt) * particle_mass; #particle/s * gram / particle = grams/s
        #https://doi.org/10.1016/j.cbpa.2020.110683
        passrate = 0.64 * 4.3 * (log10(m)) #figura 2c
        # Single day gut passage
        # passage rate of a kJ within a single particle (needs to be multiplied by kJ in stomach to get total kJ passing in a day)
         # grams/s * s/day * kJ/gram = kJ/day
        epsilon = 0.1;
        # Set up final conditions for W
        for y in 1:ymax
            if y>yc
                S[m][:, y, tmax] .= 1
            else
                S[m][:, y, tmax] .= 0
            end
        end
        
        # Main loop
        for t in (tmax-1):-1:1
            for y in (yc+1):ymax
                for x in 1:xmax
                    value = zeros(length(tid))
                    #average_gain=0.0
                    for k in 1:length(tid)
                        st = 0.0
                        for i in 1:nbins
                            gutpass = passrate * mass_dict[m][k]["gains"][i] #how many of the food is absorbed by the body
                            for j in 1:nbins
                                xp = Int64(round(x + (mass_dict[m][k]["prob"][i,j]*mass_dict[m][k]["gains"][i])-(gutpass)))
                                xp = bc(xp, 1, xmax)
                                yp = Int64(round(y + (epsilon*gutpass) - (mass_dict[m][k]["prob"][i,j]*mass_dict[m][k]["costs"][j])))
                                yp = bc(yp,yc, ymax)
                                #println(p(mass, 0.003, -1000))
                                #println(xp)
                                #println(yp)
                                st += mass_dict[m][k]["prob"][i,j]*S[m][xp, yp, t+1] 
                            end
                        end
                        value[k] = (1-d)*st
                    end
                    maxvalue = maximum(value)
                    istar = argmax(value)  
            
                    S[m][x, y, t] = maxvalue
                    #D[m][x, y, t] = tweight[istar]
                    #R[m][x, y, t] = tid[istar] 
                    T[m][x, y, t] = istar
                end
            end
        end
    end
    
    return S, T
end
println("inicio SDP")

S, T = SDP(m_min, m_max)

println("fim SDP")
# stringdataS = JSON.json(S)
# #stringdataD = JSON.json(D)
# #stringdataR = JSON.json(R)
# stringdataT = JSON.json(T)
# stringdataM = JSON.json(mass_dict)

# write the file with the stringdata variable information
open("/Users/admin/Documents/projects_yeakellab/megafauna_proj/herbivore_project/results_M_$m_min-$m_max.json", "w") do f
    JSON.print(f, mass_dict)
end

open("/Users/admin/Documents/projects_yeakellab/megafauna_proj/herbivore_project/results_S_$m_min-$m_max.json", "w") do f
    JSON.print(f, S)
end

#open("/Users/admin/Documents/Projects Yeakel Lab/Herbivores project/results_D_$m_min-$m_max.json", "w") do f
#    JSON3.write(f, stringdataD)
#end

#open("/Users/admin/Documents/Projects Yeakel Lab/Herbivores project/results_R_$m_min-$m_max.json", "w") do f
#    JSON3.write(f, stringdataR)
#end

open("/Users/admin/Documents/projects_yeakellab/megafauna_proj/herbivore_project/results_T_$m_min-$m_max.json", "w") do f
    JSON.print(f, T)
end

println("É TETRA")