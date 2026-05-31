#Source: https://ourworldindata.org/explorers/covid?tableSearch=ghana&pickerSort=asc&pickerMetric=entityName&Metric=Confirmed+deaths&Interval=Cumulative&Relative+to+population=false&country=~GHA
using DataFrames
using CSV   
using XLSX
using Dates

### Read Data 
#read_idx = 500;
limit_date = Date(2020, 11, 11); 
# Find earliest date across all dataframes
#first_date = minimum([minimum(df_D.Date), minimum(df_V.Date), minimum(df_I.Date)]); 
first_date =  Date(2020,3,12); # set to the earliest date of interest (March 12, 2020)

# Deaths (D)
#df_deaths = CSV.read("Ghana_Deaths.txt", DataFrame);#; limit=read_idx) # read only the first read_idx rows
# 1. Read the file
raw_data = read("Deaths_new.txt", String);

# 2. Extract the date block
date_match = match(r"categories:\s*\[(.*?)\]"s, raw_data);
count_match = match(r"data:\s*\[(.*?)\]"s, raw_data);

if date_match !== nothing && count_match !== nothing
    date_block = date_match.captures[1]
    count_block = count_match.captures[1]

    # 3. USE REGEX TO EXTRACT DATES
    # This looks for patterns like: Month Day, Year (e.g., Feb 15, 2020)
    # It ignores quotes, commas, and newlines entirely.
    date_elements = collect(m.match for m in eachmatch(r"[A-Z][a-z]{2} \d{1,2}, \d{4}", date_block))

    # 4. USE REGEX TO EXTRACT NUMBERS
    # This grabs any sequences of digits
    count_elements = collect(m.match for m in eachmatch(r"\d+", count_block))

    # 5. Define format and parse
    dfmt = dateformat"u d, yyyy"
    
    try
        dates = [Date(d, dfmt) for d in date_elements]
        dead_counts = parse.(Int, count_elements)

        # 6. Create DataFrame
        global df_D = DataFrame(Date = dates, Deaths = dead_counts)
    
    catch e
        # Diagnostic: If it fails, print the specific element that caused the crash
        println("Parsing failed.")
        rethrow(e)
    end
else
    println("Could not find categories or data blocks in the file.")
end
# filter out rows with dates in 2025 or 2026
df_D = filter(row ->  row.Date > first_date && row.Date <= limit_date, df_D);



# Total population (N) 
dates_N = [
    Date(2020, 12, 31),
    Date(2021, 12, 31),
    Date(2022, 12, 31)
    #=,
    Date(2023, 12, 31),
    Date(2024, 12, 31)
    =#
];


# Vaccinated (V)
df_vaccinated = CSV.read("cumulative_Vaccinated.txt", DataFrame);#; limit=read_idx) # read only the first read_idx rows
# filter out rows with dates in 2025 or 2026
df_V = filter(row ->  row.Date > first_date && row.Date <= limit_date, df_vaccinated)

# Infected (I) - Currently Infected

# 1. Read the file
raw_data = read("Currently_infected.txt", String);

# 2. Extract the date block
date_match = match(r"categories:\s*\[(.*?)\]"s, raw_data);
count_match = match(r"data:\s*\[(.*?)\]"s, raw_data);

if date_match !== nothing && count_match !== nothing
    date_block = date_match.captures[1]
    count_block = count_match.captures[1]

    # 3. USE REGEX TO EXTRACT DATES
    # This looks for patterns like: Month Day, Year (e.g., Feb 15, 2020)
    # It ignores quotes, commas, and newlines entirely.
    date_elements = collect(m.match for m in eachmatch(r"[A-Z][a-z]{2} \d{1,2}, \d{4}", date_block))

    # 4. USE REGEX TO EXTRACT NUMBERS
    # This grabs any sequences of digits
    count_elements = collect(m.match for m in eachmatch(r"\d+", count_block))

    # 5. Define format and parse
    dfmt = dateformat"u d, yyyy"
    
    try
        dates = [Date(d, dfmt) for d in date_elements]
        infected_counts = parse.(Int, count_elements)

        # 6. Create DataFrame
        global df_I = DataFrame(Date = dates, I = infected_counts)
    
    catch e
        # Diagnostic: If it fails, print the specific element that caused the crash
        println("Parsing failed.")
        rethrow(e)
    end
else
    println("Could not find categories or data blocks in the file.")
end


# Infected (I) - cumulative
#df_infected = CSV.read("cumulative_Infected.txt", DataFrame);#; limit=read_idx) # read only the first read_idx rows
# filter out rows with dates in 2025 or 2026
df_I = filter(row ->  row.Date > first_date && row.Date <= limit_date, df_I);



#################### Get Data for Reference Soltuion

# Deaths (D)
df_D.tref_D = Float64.(Dates.value.(df_D.Date .- first_date));
tref_D = df_D.tref_D;
refD = df_D.Deaths;

# Fully Vaccinated (V) - Cummulative
df_V.tref_V = Float64.(Dates.value.(df_V.Date .- first_date));
tref_V = df_V.tref_V;
refV = df_V.V_cum;

# Infected (I) - Cummulative
df_I.tref_I = Float64.(Dates.value.(df_I.Date .- first_date));
tref_I = df_I.tref_I;
refI = df_I.I;
#refI = df_I.Cases_cum;

############ # Total population (N) 

tref_N = Float64.(Dates.value.(dates_N .- first_date));
refN = Float64.([
    31887809,
    32518665,
    33149152
    #=,
    33787914,
    34427414
    =#
]);
###################

#### All models
# === PARAMETER VALUES === #
:b01,:b02,:b03,:b11,:b12,:b13,:ω
model1 = [0.0, 0.0, #k1, k2
        0.0,0.0,0.0, 6.038e-8 * 0.62811041, 6.038e-8,0.0,0.0,0.0,4.00199e-8, # =p, i.e. coefficients in aSE
          0.0,1.0,1.0,0.0,0.0, 0.0,1319.294,0.000042578,0.01000001,1.0,0.00500005,0.29,  
          0.0044,  0.0,  # no aSE
            0.0,  0.0,  0.0,  0.07142238,  0.0,  0.0,  0.1428524,  0.20005051,  0.79999398,  0.0,  
          0.01780400,  0.0805840,  0.0,  0.0,  0.92152716,  0.0,  0.41138431,  0.0,
          1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,pi]

model2 = [-1.0, 1.0, # k1, k2
    0.0, 0.0, 0.0, 1.12*0.3, 1.12*1.8, 1.12*0.3, 0.0, 0.0, 0.0, # =p, i.e. coefficients in aSE 
          1.0,1.0,1.0,1.0,0.0,  0.0001,  10000/(59*365),  1/(59*365),  0.7,  0.05,  0.86,  0.0,  
          0.018,  0.018, # no aSE
            0.4,  0.0,  0.2,  0.13,  0.0,  0.0,  0.13,  0.13978,  0.13978,  0.0,
           0.0,  0.0833,  0.0,  0.0833,  0.0,  0.0701,  0.011,  0.0, 
           1.0,1.0,1.0,pi,1.0,1.0,1.0,pi,1.0,pi]

model3 = [0.0, 0.0, #k1, k2
        0.0,0.0,0.0,0.0,0.3,0.0,0.0,0.0,0.0,  # =p, i.e. coefficients in aSE
          0.0,0.0,0.0,0.0,0.0,  0.0,  0.0,  0.0,  1.0,  0.5,  0.0,  0.0,  0.0,  0.0, # no aSE
            0.0,  0.0,  0.0,  0.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
            0.1,  0.0,  0.04,  0.0,  0.03,  0.0,  0.0, 
            1.0,1.0,1.0,pi,1.0,1.0,1.0,pi,1.0,pi]

model4 = [0.0, 0.0, #k1, k2
    0.0,0.0,0.0,0.0,2*9.1508e-9,0.0,0.0,0.0,0.0,  # =p, i.e. coefficients in aSE
          0.0,1.0,0.0,1.0,0.0,  0.0,  3251.0,  0.3349e-4,  1.0,  0.5,  0.0,  0.0,  1/16.1,  1/11.2, # no aSE
           0.0,  0.0,  0.0,  1/4.2,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  2/30,  0.0,  2/1.5,  0.0,  1/11.5,  0.0426,  0.0, 
          1.0,1.0,1.0,pi,1.0,1.0,1.0,pi,1.0,pi]

model5 = [1.0, 0.0, #k1, k2
    0.0,0.0,0.999, 0.0, 0.999, 0.0,0.0,0.0,0.0,  # =p, i.e. coefficients in aSE
          1.0,1.0,0.0,0.0,1.0,  1/40,  0.46,  0.0991,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0, # no aSE
          0.4,  0.0,  0.3002,  0.001,  0.0,  0.2,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.005,  0.0,  0.0,  0.0,  0.0,  1/14,
          1.0,1.0,1.0,pi,1.0,1.0,1.0,pi,1.0,pi]
all_models = [model1, model2, model3, model4, model5];


#=
using Plots
# Switch the backend to PyPlot
pyplot()

# 1. Plot Cumulative Deaths
plot_deaths = plot(
    tref_D, 
    refD, 
    label="Cumulative Deaths", 
    color=:crimson, 
    linewidth=2.5, 
    marker=:circle, 
    markersize=3,
    xlabel="Days since March 12, 2020", 
    ylabel="Number of People",
    title="Ghana COVID-19: Cumulative Deaths",
    grid=:true,
    legend=:topleft,
    left_margin=15Plots.mm # Adds explicit margin padding so y-label isn't cut off
)

# 2. Plot Currently Infected
plot_infected = plot(
    tref_I, 
    refI, 
    label="Active Cases", 
    color=:royalblue, 
    linewidth=2.5, 
    marker=:square, 
    markersize=3,
    xlabel="Days since March 12, 2020", 
    ylabel="Number of People",
    title="Ghana COVID-19: Active Infections",
    grid=:true,
    legend=:topright,
    left_margin=15Plots.mm # Adds padding to prevent the y-axis text collision
)

# 3. Combine with a tighter layout adjustment
combined_plot = plot(plot_deaths, plot_infected, layout=(1, 2), size=(1100, 450))

display(combined_plot)
=#