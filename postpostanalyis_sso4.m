function postpostanalyis_sso4

%nested cell arrays: array{n}{m}
%find non-empty elements: test3=find(~cellfun(@isempty,test));

%% Input variables

filelist={'S01-processed(80-40)-sso4.mat','S02-processed(80-40)-sso4.mat','S03-processed(80-40)-sso4.mat','S04-processed(80-40)-sso4.mat','S05-processed(80-40)-sso4.mat','S06-processed(80-40)-sso4.mat'};
SequencesToTest=10; % How many sequences should be traced through the selection rounds
SimilaritySeqs=50; % How many shortened core sequences should be compared with everything for mutations
AlignmentThreshold=0.95; % Cutoff for whether sequence counts as mutated version or not.

%% Data loading

for(a=1:1:length(filelist))
    load(filelist{a});
    SeqCores{a}=SortedSeqCore;
    Seqs{a}=SortedSeq;
    Percentages{a}=EnrichmentPercentage;
    Counts{a}=EnrichmentCount;
    
    varstoclear={'SortedSeqCore','SortedSeq','EnrichmentPercentage','EnrichmentCount'};
    clear(varstoclear{:});
end

%% Counting occurences in each round 

for(a=1:1:length(filelist))

    for(b=1:1:length(filelist))
        
        for(c=1:1:SequencesToTest)
            
            %checks here that c is within bounds.
            if(c<=length(Seqs{a}))
               
                temp=strfind(Seqs{b},Seqs{a}{c});             
                positiontemp=find(~cellfun(@isempty,temp));
            
            else
                
                temp=[];
                positiontemp=[];
                
            end
            
            
            if(isempty(positiontemp))
                result(b,c)=0;
            else
                result(b,c)=positiontemp(1,1);
            end
            
        end
    end

    AllResults{a}=result;
    
end

%% Re-Sorting shorted sequences

tic

for(a=1:1:length(filelist))
    
    [CoreSeq2,CorePercentage2]=EliminateDuplicates(SeqCores{a},Percentages{a});

        
    [temp, sortedindex]=sort(CorePercentage2,'descend');
    
    CorePercentage{a}=CorePercentage2(sortedindex);
    SeqCores{a}=CoreSeq2(sortedindex);
        
end

toc

%% Finding related sequences
% align sequence with itself, then use set percentage of score as cutoff



for(a=1:1:SimilaritySeqs)
    
    TempStruct=localalign(SeqCores{end}{a},SeqCores{end}{a});

    Threshold=TempStruct.Score*AlignmentThreshold;
    
    counter=1;
    
    SimilaritySequences{a}{counter}=SeqCores{end}(a);
    SimilarityPercentages{a}(counter)=CorePercentage{end}(a);
        
    for(b=a+1:1:length(SeqCores{end}))
       
        if(~isempty(SeqCores{end}{b}))
        
            TempStruct=localalign(SeqCores{end}{a},SeqCores{end}{b});
            
            if(TempStruct.Score>=Threshold)
                
                counter=counter+1;
                
                SimilaritySequences{a}{counter}=SeqCores{end}(b);
                SimilarityPercentages{a}(counter)=CorePercentage{end}(b);
            end
        
        end
        
    end
    
end

%% Saving

for(a=1:1:length(filelist))
    Filename=['TopSequencesofRound' num2str(a) 'throughoutSelectionTemp'];
        
    [filehandle, message]=fopen([Filename '.txt'], 'w');

    if(filehandle==-1)
        error(['Error opening file ' Filename ': ' message]);
    end
        
    % write header
    fprintf(filehandle, 'Sequence');
    
    for(b=1:1:length(filelist))
        fprintf(filehandle, '\tRound%u',b);
    end
    
    fprintf(filehandle, '\n');
    
    for(b=1:1:SequencesToTest)
        fprintf(filehandle,'%s', Seqs{a}{b});
        
        for(c=1:1:length(filelist))
            temp=AllResults{a}(c,b);
            
            if(temp~=0)
                returnvalue=Percentages{c}(temp);
            else
                returnvalue=0;
            end
            
            
            fprintf(filehandle, '\t%g',returnvalue);
        end
        
        fprintf(filehandle, '\n');
        
    end
    
    fclose(filehandle);

    
end

for(a=1:1:SimilaritySeqs)

Filename=['SimilarToSequence' num2str(a) 'AlignmentThreshold' num2str(AlignmentThreshold)];
        
    [filehandle, message]=fopen([Filename '.txt'], 'w');

    if(filehandle==-1)
        error(['Error opening file ' Filename ': ' message]);
    end
        
    % write header
    fprintf(filehandle, 'Sequence\tFrequency(Percent)\n');

     
    for(b=1:1:length(SimilaritySequences{a}))
                
        fprintf(filehandle, '%s\t%g\n', SimilaritySequences{a}{b}{:},SimilarityPercentages{a}(b) );

    end
    
    fclose(filehandle);

    
end
        
