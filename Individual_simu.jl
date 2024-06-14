using Distributions
using Random
using JSON
using LinearAlgebra

# configurations=10

# function read_json(results_S.json)
#     open(stringdataS,"r") do f
#         return JSON.parse(f)
#     end
# end
function mean_retention_time(mass)
    mean_retention_time = (30.3 * (mass)^0.109)
    return mean_retention_time * 60 * 60 # [hr -> s]

end

function gut_volume_g(mass)
    capacity = (0.030 * (mass)^0.881);
    return capacity  #[kg -> g]
    #return capacity
end



function mean_particle_mass(mass)
    # from "comparative chewing efficiency"
    # mean particle size in mm
    # mass in g

    mass_g = mass * 1000; # convert mass kg -> g
    #mean_particle_size = 0.0;
    mean_particle_size = (6.61 * (mass_g)^0.26)
    #mean_particle_size = (6.61 * (mass)^0.26)
    volume = (4/3) * pi * (1/2 * mean_particle_size)^3; # [mm^3]
    particle_mass = 0.0004 * volume; # [g/mm^3 * mm^3 = g]

    return particle_mass

end

function bc(yvalue, yc, ymax)
    return clamp(yvalue, yc, ymax)
end ####### FUNÇÃO DE CORTE #############



S = Dict()
S = JSON.parsefile("/Users/admin/Documents/Projects Yeakel Lab/Herbivores project/results_S_100-110.json")
#     global S
#     # dicttxt = readall(f)  # file information to string
#     S=JSON.parse(f)  # parse and transform data
# end

# D = Dict()
# open("/Users/admin/Documents/Projects Yeakel Lab/Herbivores project/results_D_30-32.json", "r") do f
#     global D
#     # dicttxt = readall(f)  # file information to string
#     D=JSON.parse(f)  # parse and transform data
# end

# R = Dict()
# open("/Users/admin/Documents/Projects Yeakel Lab/Herbivores project/results_R_30-32.json", "r") do f
#     global R
#     # dicttxt = readall(f)  # file information to string
#     R=JSON.parse(f)  # parse and transform data
# end

T = Dict()
T = JSON.parsefile("/Users/admin/Documents/Projects Yeakel Lab/Herbivores project/results_T_100-110.json")
#     global T
#     dicttxt = read(f, String)  # file information to string
#     T=JSON.parse(dicttxt)  # parse and transform data
# end

#println(typeof(T))

mass_dict = Dict()
open("/Users/admin/Documents/Projects Yeakel Lab/Herbivores project/results_M_100-110.json", "r") do f 
     global mass_dict
     #dicttxt = read(f)  # file information to string
     mass_dict=JSON.parse(f)  # parse and transform data
 end

#println(typeof(mass_dict))
# mass_dict = Dict()
# open("results_massdict_10-30.json", "r") do f
#     global mass_dict
#     # dicttxt = readall(f)  # file information to string
#     mass_dict=JSON.parse(f)  # parse and transform data
# end
tmax=29
configurations=100
dead = 0
final_state = Dict()
for m in 100:110
    final_state[m]= Dict()
    ymax = Int64(round(4.87* 0.02 * (m^1.19)))#max fat storage in Mcal
    yc = Int64(round(4.87*0.3*0.02 * (m^1.19)))
    xmax = Int64(round(gut_volume_g(m)*4.87));
    for k in 1:configurations
    #for t in 1:length(S[m])
        for t in 1:tmax
            if t==1
                initial_cond_y= Int64(round(rand(Uniform(yc+1, ymax))))
                initial_cond_x= Int64(round(rand(Uniform(1, xmax))))
            end
            # final_state[m][t]= zeros(Float64, tmax, configurations)
            #println(initial_cond_x)
            #println(initial_cond_y)
            decision = Int64(T["$m"][t][initial_cond_y][initial_cond_x])
            #println(t, decision)
            gains = zeros(Float64,50)#nbins=50
            r= rand(Uniform(0,1))
            #println(r)
            if r>S["$m"][t][initial_cond_y][initial_cond_x]
                #final_state[m][t]= zeros(Float64, configurations)
                #final_state[m][t][k]=decision
                #println(r, S["$m"][t][initial_cond_y][initial_cond_x])
                global dead+=1
                #println("no céu tem pão????")
                break
            
            else
                for i in 1:length(mass_dict["$m"]["$decision"]["prob"])
                    gains[i] = sum(mass_dict["$m"]["$decision"]["prob"][i])#soma de uma linha
                end

                gains = gains./sum(gains) #normalizar a probabilidade (probabilidade de ganho)
                #println("gains =", sum(gains))
                cleber = Multinomial(1, gains)
                jequiti = rand(cleber)
                jequiti_index= findfirst(iten->iten==1, jequiti)
                gain = mass_dict["$m"]["$decision"]["gains"][jequiti_index]
                #println("olha o roda-roda jequiti=",gain)
                cost = mass_dict["$m"]["$decision"]["prob"][jequiti_index]
                costs = cost./sum(cost) 
                cleber_c =  Multinomial(1, costs)
                gugu = rand(cleber_c)
                gugu_index= findfirst(iten->iten==1, gugu)
                cost =  mass_dict["$m"]["$decision"]["costs"][gugu_index]
                #println("o show do gugu!!! =", cost)
                ########################################## new y and x ########################################################################
                #mrt = mean_retention_time(m); # seconds / PARTICLE
                
                #particle_mass = mean_particle_mass(m); #gram / particle
                
                # Passage rate of food (rate of flow from gut to body)
                passrate = 0.64 * 4.3 * (log10(m)); #particle/s * gram / particle = grams/s
                
                # Single day gut passage
                # passage rate of a kJ within a single particle (needs to be multiplied by kJ in stomach to get total kJ passing in a day)
                # grams/s * s/day * kJ/gram = kJ/day
                epsilon = 0.1;
                gutpass = passrate * gain #how many of the food is absorbed by the body
                
                new_xp = Int64(round(initial_cond_x + gain -(gutpass)))
                new_xp = bc(new_xp, 1, xmax)
                new_yp = Int64(round(initial_cond_y  + (epsilon*gutpass) - cost))
                new_yp = bc(new_yp,yc, ymax)
                global initial_cond_x = new_xp
                global initial_cond_y = new_yp

                if t==10
                    final_state[m][t]= zeros(Float64, configurations)
                    final_state[m][t][k]=decision #(?)
                elseif t==20
                    final_state[m][t]= zeros(Float64, configurations)
                    final_state[m][t][k]=decision #(?)
                elseif t==29
                    final_state[m][t]= zeros(Float64, configurations)
                    final_state[m][t][k]=decision #(?)
                end
            end
        end
    end    
end 
println(dead)
