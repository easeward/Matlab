%This is the first script that should be run and generates the
%distributions for the paramters you input.Best not to run it much bigger
%than 10x10x10x10. 
%The outputs that are most important are total_cell_T_array_nucleus and 
% total_cell_T_array_cytoplasm as these are the distributions in the
% nucleus and cytoplasm that will be compared to real data. 
tic %Starts timing how long the whole run takes
% DON'T FORGET TO CHANGE THE NAME OF THIS RUN AT
% THE END OF THE SCRIPT TO SOMETHING LOGICAL AND INFORMATIVE

%Use the automated parameter filler below if you want evenly spaced inputs
% PARAMETER FILLER
% from1 = 0.1;
% interval1 = 6;
% to1 = 40;
% from2 = 0.1;
% interval2 = 8;
% to2 = 73;
% from3 = 0.01;
% interval3 = 0.5;
% to3 = 5;
% from4 = 0;
% interval4 = 6;
% to4 = 30;
%av_time_off_input = [from1:interval1:to1];
%av_time_on_input = [from2:interval2:to2];
%c_input = [from3:interval3:to3];
%elongation_input = [from4:interval4:to4];

%Use this to manually input your parameters
av_time_off_input = [ 0.1 0.2 0.4 0.8 1 2 5 10 20 40];
av_time_on_input = [ 0.1 0.2 0.4 0.8 1 2 5 10 20 40 60 85];
c_input = [ 0.01 0.02 0.05 0.1 0.2 0.5 1 2 4 6];
elongation_input = [0.005 0.09 0.11 0.13 0.17 0.19 0.3 0.5 1.0 2.0 4.0];
%can add another inupt
%Don't forget that a = 1/average time off so increasingly big average time
%off values will make very little difference to the value of a. Likewise b.

iteration_av_time_off = 1; %This starts the iteration counter.
max_transcripts_all_nucleus = 0; %Before model starts there are no transcripts
max_transcripts_all_cytoplasm = 0;
max_transcripts_nucleus = 19; %Cuts off data at a biologically relevant size. Change if you want bigger or smaller
max_transcripts_cytoplasm = 29; %Likewise.

%THIS IS FOR TESTING 1 SET OF INPUTS AT A TIME. 1X1X1X1 MATRIX. COMMENT OUT
%FOR REAL RUN OF MODEL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% av_time_off_input = 1.4; %average time that the gene is off (mins)
% av_time_on_input = 8.7; %average time the gene is on (mins)
% c_input = 0.2; %initiation rate mins^-1. if c/b > 1 there is bursting.
%  length of the gene in kb
% elongation_input = 2;

for off = av_time_off_input; 
    %iteration_av_time_off; %Take away the ; if you want to keep track of this loop
    iteration_av_time_on = 1; %Resets to 1 after each loop is completed
    for on = av_time_on_input;
    %iteration_av_time_on ;
    iteration_c_values = 1; %Resets to 1 after each loop is completed
        for c = c_input;
            %iteration_c_values;
            iteration_elongation = 1;%Resets to 1 after each loop is completed
            for elongation = elongation_input;
          
            iteration_av_time_off %Let's you keep track fo model progression
            iteration_av_time_on
            iteration_c_values
            iteration_elongation
  % parameters = [av_time_off av_time_on c gene_length elongation_rate generation_time number_of_cells]
            parameters = [off on c 1.8 elongation 85 1000]; 
            %must have at least 2 cells (final value) if you want to be
            %able to average distributions.!!!!!!!!!
            %Don't forget to change the gene length for WT (1.8kb) and
            %mutants (0.803kp)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 %EXAMPLES OF INPUTS 
%av_time_off = 1.4; %average time that the gene is off (mins)
%av_time_on = 8.7; %average time the gene is on (mins)
%c = 0.2; %initiation rate mins^-1. if c/b > 1 there is bursting.
%gene_length = 1.8; %length of the gene in kb
%elongation_rate = 2; %elongation rate in kb min^-1
number_of_cells = parameters(7);
for z = 1:number_of_cells %number of cells
%THE MODEL
%INITIAL CONDITIONS (INITIALISATION STEP)
    av_time_off= parameters(1);
    av_time_on = parameters(2);
    c = parameters(3);
    gene_length = parameters(4);
    elongation_rate = parameters(5);
    generation_time = parameters(6);
    T_nucleus = 0; %no transcripts to start with
    T_cytoplasm = 0; %no transcripts to start with
    a  = (1/av_time_off); %rate of activation min^-1 
    b =(1/av_time_on); %rate of inactivation min^-1 
    state = 0; % state = 0 means inactive state = 1 means active gene.Always start inactive.
    nucleus_removal_time = (gene_length/elongation_rate); %time (mins) before transcript leaves nucleus 
    on_state_rates = [b c]; %chance of b, inactivation, or c, initiation
    t_final = generation_time; %Model stops after generation time reached (85 mins) can vary if necessary 
    

% EVENTS AND TIME BETWEEN EVENTS (MONTE CARLO STEP)
    
    time_increase_nucleus = []; %empty vector times (mins) between events
    initiation = []; %empty vector marking initation of transcription events
    current_time_nucleus = 0;
    i = 0;
    finaltime = -1;
    while finaltime~=i
    if current_time_nucleus <= t_final;
        i = i+1;
        if state == 0; %If the gene is inactive the only scenario is activation
            %k = race_conditions(a); %Add in race script if more than one
            %scenario is wanted at this stage.
            tau = log(1/rand)*(1/a); %time step tau(mins) until activation
            time_increase_nucleus(i) = tau;
            state = 1; %Gene is now active
            initiation(i) = 0; %gene becomes active but no initiation yet.
        elseif state ==1;
            tau = log(1/rand)*(1/(sum(on_state_rates))); %Time step until next event
            time_increase_nucleus(i) = tau;
            %"RACE" THE TWO SCENARIOS TO FIND NEXT EVENT
            %on_state_rates is a vector that contains probabilities of all possible
            %events eg. [0.1 0.2];
            cumulativeProbs = cumsum(on_state_rates);%vector of cumulative probs of b or c eg. [0.1 0.3]
            total_on_state_rates = sum(on_state_rates); %eg. 0.3
            U = total_on_state_rates*rand; %rand is a random number 0-1 with normal distribution.
            %multiplying by total_on_state_rates means the random numbers
            %will be generated in that range.
            sample = 1; %start with the assumption that the initial scenario (b in this case) will be chosen. 
            while U > cumulativeProbs(sample)
                sample = sample + 1; %If the random number is greater than the initial scenario probability, the second scenario is chosen. 
            end
            
            x = sample; %eg 1 or 2 if b or c chosen.
            if x == 1;
                %disp('return to inactive')
                T_nucleus = T_nucleus; %The number of transcripts remains the same.
                state = 0; %Gene returns to inactive
                initiation(i) = 0; %No initiation of transcription
           
            elseif x ==2;
                %disp('return to choice state & initiaion set off')
                T_nucleus = T_nucleus + 1; % adds transcript number to T_nucleus.
                state = 1; %gene remains active
                initiation(i) = 1; %Transcription is initiated
            else
                disp('modify code to deal with more than 2 probabilites')
            end
        else
            disp ('error')
        end
        time = cumsum(time_increase_nucleus); %Keeps track of time
        current_time_nucleus = time(end); %current time is last time recorded
    else 
            i = finaltime; %stop running the model when generation time reached. 
    end
    
    end
    %TAKE TIMESTEPS AND TURN INTO TOTAL TIME SPENT WITH THAT NUMBER OF
    %TRANSCRIPT

    %making vector with time of transcript initiation or decay on top
    %and total level of transcript at that time on the bottom. 
    event_times_nucleus = cumsum(time_increase_nucleus); %times of any event happening
    %initiation; %vector of times when initiation starts
    transcripts_cumulative = cumsum(initiation); %total number of transcripts as time increases
    initiation_times = event_times_nucleus(logical(initiation)); %times of initiations only
    removal_times = initiation_times + nucleus_removal_time; %times of decay set at d after initiation
    removal_times_real = removal_times(removal_times <= generation_time); %only choses those within the time frame
    cytoplasm_entry_times = [removal_times_real, generation_time]; %add 85 to end to cap
    Transcript_Half_life = 19; %minutes Perhaps loop this
    decay_rate = log(0.5)/Transcript_Half_life; %minutes^-1
    transcript_destruction = -(decay_rate); %to make it a positive number
    
    
    %MODEL DECAY OF TRANSCRIPTS IN THE CYTOPLASM
    time_increase_cytoplasm = []; %empty vector times (mins) between events
    transcript_destruction_times = []; %empty vector marking initation events
    cytoplasm_transcripts = [];
    current_time_cytoplasm = 0;
    time_increase_cytoplasm(1) = current_time_cytoplasm;
    cytoplasm_transcripts(1) = 0; % change this if you want transcripts at time zero
    %decisions_cytoplasm = [(transcript_destruction)]; %at the moment only
    %decay, if more options add here.
    current_time_cytoplasm = cytoplasm_entry_times(1);
    time_increase_cytoplasm(2) = current_time_cytoplasm;
    cytoplasm_transcripts(2) = 1;
    i = 0;
    decay_step = 2; %the cytoplasm_entry time we are comparing against
    while j~=i
        if current_time_cytoplasm < t_final;
            i = i+1;
            %cytoplasm_entry_times(decay_step);
            %cytoplasm_transcripts(end);
            tau = log(1/rand)*(1/((transcript_destruction)*(cytoplasm_transcripts(end)))); %time step tau(mins) until transcript is removed (more likely with more transcripts)
            if tau >=generation_time
                tau = generation_time; %stops the tau going to infinity
            end
            %(cytoplasm_entry_times(decay_step) - current_time_cytoplasm);
                if tau <= (cytoplasm_entry_times(decay_step) - current_time_cytoplasm);
                    %If a transcript will decay before a new one enters do
                    %this...
                    current_time_cytoplasm = current_time_cytoplasm + tau; %Update current time
                    time_increase_cytoplasm(i+2) = current_time_cytoplasm;
                    %time_increase_cytoplasm;
                    transcript_destruction_times(i) = current_time_cytoplasm;%Record time when decayed
                    cytoplasm_transcripts(i+2) = (cytoplasm_transcripts(end) -1); %Remove 1 transcript
                    %cytoplasm_transcripts;
                    %current_time_cytoplasm;
                elseif tau >= (cytoplasm_entry_times(decay_step) - current_time_cytoplasm);
                    %if a new transcript will enter before the next decay
                    %event do this....
                    current_time_cytoplasm = cytoplasm_entry_times(decay_step); %Update time to next entry time
                    time_increase_cytoplasm(i+2) = current_time_cytoplasm;
                    %time_increase_cytoplasm;
                    if current_time_cytoplasm == generation_time;
                        %If the cell has reached the generation cap keep the
                    %number of transcripts constant to add end point
                        cytoplasm_transcripts(i+2) = (cytoplasm_transcripts(end));
                    else
                    cytoplasm_transcripts(i+2) = (cytoplasm_transcripts(end) +1); %increase number of transcripts by 1
                    %cytoplasm_transcripts;
                    decay_step = decay_step +1;
                    %current_time_cytoplasm;
                    end
                else
                    disp 'Something is amiss'
                end
        else 
            i = j;
        end
    
    end
    cytoplasm_transcripts; %Number of transcripts in cytoplasm
    time_increase_cytoplasm; %Timing of events in cytoplasm
    
    %SORT DATA TO GIVE VECTOR OF [NUMBER OF TRANSCRIPTS, TIME POINTS]
    event_transcript_changes_unsorted = [0 initiation_times removal_times_real];
    T_decayed = length(removal_times_real);%adds time point 0
    T_changing = zeros(1,(T_nucleus+T_decayed+1)); %empty vector to fill
    plus = length(initiation_times);
    T_changing(1) = 0; %makes first time point no change
    T_changing(2:(plus+1)) = 1; %Add a transcript for every initiation event
    T_changing((plus+2) : length(event_transcript_changes_unsorted)) = -1; %Remove transcript for every removal event
    time_over_event = [event_transcript_changes_unsorted ; T_changing];
    sorted = sortrows(time_over_event',1)';%sorts changes according to time

    %CAP THE READINGS AT GENERATION TIME 85 minutes
    generation_cap_is = [sorted(1, :) <= generation_time];
    generation_cap_events = generation_cap_is(generation_cap_is ==1);
    cap_length = length(generation_cap_events);
    generation_cap_sorted = sorted(1:2, 1:cap_length);
    ending = [generation_time; 0];
    generation_cap_sorted_good = [generation_cap_sorted ending];


    T_changing_times = generation_cap_sorted_good(1, 1:end); %
    T_cumulative_nucleus = cumsum(generation_cap_sorted_good(2, 1:end)); %
    data_nucleus = [T_changing_times; T_cumulative_nucleus];
    
    %FIND TIME SPENT WITH EACH NUMBER OF TRANSCRIPTS (NUCLEUS)
    iterations_nuclues = length(T_changing_times);
    tau_event = zeros(1,(iterations_nuclues - 1));
    for i = 1:(iterations_nuclues - 1)
        tau_event (i) = generation_cap_sorted_good(1, i+1) - generation_cap_sorted_good(1,i);
    end
    tau_events_nucleus = [tau_event 0];% add 0 so it reads that it is at 85 for 0 time.

    times_T = [ T_cumulative_nucleus; tau_events_nucleus];
    
    %FIND TIME SPENT WITH EACH NUMBER OF TRANSCRIPTS (CYTOPLASM)
        iterations_cytoplasm = length(time_increase_cytoplasm);
        tau_event_cytoplasm = zeros(1,(iterations_cytoplasm - 1));
    for i = 1:(iterations_cytoplasm - 1)
        tau_event_cytoplasm (i) = time_increase_cytoplasm(i+1) - time_increase_cytoplasm(i);
    end
    tau_events_cytoplasm = [tau_event_cytoplasm 0];% add 0 so it reads that it is at 85 for 0 time.

    times_T_cytoplasm = [ cytoplasm_transcripts; tau_events_cytoplasm];
    

    %most of the tau_event times are the rate of export. Not a bug

    ordered_nucleus = sortrows(times_T',1)';
    ordered_cytoplasm = sortrows(times_T_cytoplasm',1)';
    total_time_at_T_nucleus = zeros(1, (max_transcripts_nucleus+1));%empty matrix to fill
    total_time_at_T_cytoplasm = zeros(1, (max_transcripts_cytoplasm+1));

    for i = 0:(max_transcripts_nucleus) %This is where I cap it to biologically relevant transcript numbers
        x =  ordered_nucleus(1,:)==i;
        x2 = zeros(1,length(x));
        x3 = [x2;x];
        i_total_time = ordered_nucleus(logical(x3));
        y = sum(i_total_time);
        total_time_at_T_nucleus(i+1) = y;
    
    end
    
    total_time_at_T_nucleus %We love this output. It is what we have been searching for.
    total_time_at_Ta_nucleus{z}= [total_time_at_T_nucleus]; %This keeps track of it for each cell

    
    for i = 0:(max_transcripts_cytoplasm) %This is where I cap it to be biologically relevant transcript numbers
        x =  ordered_cytoplasm(1,:)==i;
        x2 = zeros(1,length(x));
        x3 = [x2;x];
        i_total_time_cytoplasm = ordered_cytoplasm(logical(x3));
        y = sum(i_total_time_cytoplasm);
        total_time_at_T_cytoplasm(i+1) = y;
    
    end
    
    total_time_at_T_cytoplasm %this is also a great output. 
    total_time_at_Ta_cytoplasm{z}= [total_time_at_T_cytoplasm]; %This keeps track of it for each cell

z = z+1;
end

 T_number_array_nucleus = [0 : max_transcripts_nucleus];
 cell_T_array_nucleus = zeros(number_of_cells, length(T_number_array_nucleus));
 T_number_array_cytoplasm = [0: max_transcripts_cytoplasm];
 cell_T_array_cytoplasm = zeros(number_of_cells, length(T_number_array_cytoplasm));
    
%total_time_at_Ta{1}

for z=1:number_of_cells
    length_times = length(total_time_at_Ta_nucleus{z});
    length_array = length(T_number_array_nucleus);
    length_difference = length_array-length_times;
       if length_times == length_array
            cell_T_array_nucleus(z, 1:end) = [total_time_at_Ta_nucleus{z}];
    
       else 
           cell_T_array_nucleus(z, 1:end) = [total_time_at_Ta_nucleus{z}, zeros(1,length_difference)];
      
           
       end
end
%now for the cytoplasm
for z=1:number_of_cells
    length_times_cytoplasm = length(total_time_at_Ta_cytoplasm{z});
    length_array_cytoplasm = length(T_number_array_cytoplasm);
    length_difference_cytoplasm = length_array_cytoplasm-length_times_cytoplasm;
       if length_times_cytoplasm == length_array_cytoplasm
            cell_T_array_cytoplasm(z, 1:end) = [total_time_at_Ta_cytoplasm{z}];
    
       else 
           cell_T_array_cytoplasm(z, 1:end) = [total_time_at_Ta_cytoplasm{z}, zeros(1,length_difference_cytoplasm)];
      
           
       end
end

add_nucleus = sum(cell_T_array_nucleus); %Adds the values for all of the cells at that time (unually 1000 cells)
add_cytoplasm = sum(cell_T_array_cytoplasm);
average_nucleus = add_nucleus/number_of_cells %Averages the distributions
%transcript_number_nucleus = [0: max(max_transcripts_all_nucleus)];
average_cytoplasm = add_cytoplasm/number_of_cells
%transcript_number_cytoplasm = [0: max(max_transcripts_all_cytoplasm)];
%normal = average_nucleus/generation_time;

total_cell_T_array_nucleus{iteration_av_time_off, iteration_av_time_on, iteration_c_values,iteration_elongation} = average_nucleus;

total_cell_T_array_cytoplasm{iteration_av_time_off, iteration_av_time_on, iteration_c_values, iteration_elongation} = average_cytoplasm;

                iteration_elongation = iteration_elongation +1;
            end

            iteration_c_values = iteration_c_values+1;
            end
    iteration_av_time_on = iteration_av_time_on+1;
    end
    iteration_av_time_off = iteration_av_time_off+1;
end

toc
save('INSERTNAMEYOUWANT.mat') %THIS SAVES YOUR MATRIX OF DISTRIBUTIONS AND IS VERY IMPORTANT
%CHANGE IT TO SOMETHING NEW AND INFORMATIVE EACH TIME YOU RUN THE SCRIPT

%If you are developing this script and want to look at the non-capped
%version where transcripts are recorded even if they go up to >200
%transcripts then look at the script 'nucleus_cytoplasm_together.m' or
%email me (emily.seward@biodtp.ox.ac.uk). 