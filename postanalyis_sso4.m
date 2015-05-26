function postanalyis_sso4
load ('S01_raw.mat');
%clear('Sequences');

OutputFilename='S01-processed(80-40)-sso4';

%How Many Sequences to Analyze
TopSequencesToAnalyze=50000;

%Status for each sequence and static variables
Status=zeros(TopSequencesToAnalyze,1);
OK=0;
WrongRV=1;
WrongFW=2;
Merged=3;

ForwardPrimerSeq='GGACGACCTAAGGCAAACGCTATGGTCGTTAGTATGGTCGTTA' %106 max cutoff
ReversePrimerSeq='CCAGTCTCAACGTCGAGTTACGAAGA' %68 max cutoff
FW_CutOff=80 %Alignment score to exclude bad forward alignment
RV_CutOff=40 %Alignment score to exclude bad reverse alignment

%read top sequences
for(a=1:1:TopSequencesToAnalyze)
    Sequences{a}=SequencesSorted{order(a)};
    Freq(a)=frequencies(order(a));
    Count(a)=SequenceCount(order(a));
    OriginalSequences{a}=SequencesSorted{order(a)};
    OriginalFreq(a)=frequencies(order(a));
    OriginalCount(a)=SequenceCount(order(a));
end

%trim after reverse primer, sort out seq. w/ wrong reverse primer

for(a=1:1:TopSequencesToAnalyze)
    TempStruct=localalign(Sequences{a},ReversePrimerSeq);
    if(TempStruct.Score<RV_CutOff)
        Status(a)=WrongRV;
        continue;
    end
    temp=Sequences{a};
    Sequences{a}=temp(1:TempStruct.Stop(1));
end
    
%sort out sequences with wrong forward primer

for(a=1:1:TopSequencesToAnalyze)
    TempStruct=localalign(Sequences{a},ForwardPrimerSeq);
    if(TempStruct.Score<FW_CutOff)
        Status(a)=WrongFW;
        continue;
    end
    temp=Sequences{a};
    Sequences{a}=temp(TempStruct.Start(1):end);
end


%Status Output

TotalWrongFw=sum(Status==WrongFW)
TotalWrongRv=sum(Status==WrongRV)

%Merge Sequences, loop should be fast enough up to ~10000 Sequences
tic;
for(a=1:1:TopSequencesToAnalyze)
    
    if(Status(a)==OK)
      
        for(b=a+1:1:TopSequencesToAnalyze)
            if(strcmp(Sequences{a},Sequences{b}))
                %Sequences{b}=[];
                Status(b)=Merged;
                Freq(a)=Freq(a)+Freq(b);
                Count(a)=Count(a)+Count(b);
            end
            
        end
    end
end
toc

%some statistics on filtering
TotalFilteredNumber=sum(Status==OK)
TotalFilteredFreq=sum(Freq(Status==OK))
TotalMergedFreq=sum(Freq(Status==Merged))
TotalDiscardedFreq=sum(Freq(Status==WrongFW | Status==WrongRV))

%seperate filtered sequences
FilteredSequences=Sequences(Status==OK);
FilteredSequencesCount=Count(Status==OK);
FilteredSequencesFreq=Freq(Status==OK);
DiscardedSequences=OriginalSequences(Status==WrongFW | Status==WrongRV);
DiscardedSequencesCount=OriginalCount(Status==WrongFW | Status==WrongRV);
DiscardedSequencesFreq=OriginalFreq(Status==WrongFW | Status==WrongRV);




%Motif finding
RandomRegionLimits=zeros(TotalFilteredNumber,2);
%alignment is fast, redo it instead of messing with previous loop
%define random region
for(a=1:1:TotalFilteredNumber)
    TempStruct=localalign(FilteredSequences{a},ForwardPrimerSeq);
    RandomRegionLimits(a,1)=TempStruct.Stop(1)+1;
    TempStruct=localalign(FilteredSequences{a},ReversePrimerSeq);
    RandomRegionLimits(a,2)=TempStruct.Start(1)-1;
end

%these are the three motifs and respective cutoff values
Motif1A='CGATC'; %Max: 14.6667, 1 mismatch=12
Motif1ACutOff=11.9;
Motif1B='GATCG'; %Max: 13, 1 mismatch=11.3333
Motif1BCutOff=11.2;
Motif2='TGGAAACA'; %Max: 18, 1 mismatch=16.3333
Motif2CutOff=16.2;
Motif3='ATGTT'; %max: 9.33, 1 mismatch=7 to 7.666
Motif3CutOff=7;
Motif4='TC';
MotifPresent=zeros(TotalFilteredNumber,5);
LibraryClass=zeros(TotalFilteredNumber,1);


for(a=1:1:TotalFilteredNumber)
    
    %check for the three motifs
    temp=FilteredSequences{a};
    
    if(RandomRegionLimits(a,1)+15<=length(temp))
        TempStruct=localalign(temp(RandomRegionLimits(a,1):RandomRegionLimits(a,1)+15),Motif1A);
        if(TempStruct.Score>=Motif1ACutOff)
            MotifPresent(a,1)=1;
        end
    end
    
    if(RandomRegionLimits(a,2)-15>=1)
        TempStruct=localalign(temp(RandomRegionLimits(a,2)-15:RandomRegionLimits(a,2)),Motif1B);
        if(TempStruct.Score>=Motif1BCutOff)
            MotifPresent(a,2)=1;
        end
    end
    
    if(RandomRegionLimits(a,1)<RandomRegionLimits(a,2))
        TempStruct=localalign(temp(RandomRegionLimits(a,1):RandomRegionLimits(a,2)),Motif2);
        if(TempStruct.Score>=Motif2CutOff)
            MotifPresent(a,3)=1;
        end
    end
    
    if(RandomRegionLimits(a,1)<RandomRegionLimits(a,2))
        TempStruct=localalign(temp(RandomRegionLimits(a,1):end),Motif3);
        if(TempStruct.Score>=Motif3CutOff)
            MotifPresent(a,4)=1;
            
            if(TempStruct.Stop(1)+10<RandomRegionLimits(a,2) && ~isempty(strfind(temp(TempStruct.Stop(1)+10+1:RandomRegionLimits(a,2)), Motif4)))
                MotifPresent(a,5)=1;
            end
        end
    end
        

    %group into library classes according to mofifs present
    if(MotifPresent(a,1) && MotifPresent(a,2) && MotifPresent(a,3))
        LibraryClass(a)=1;
        continue;
    end
    
    if(MotifPresent(a,1) && MotifPresent(a,2))
        LibraryClass(a)=2;
        continue;
    end
    
    LibraryClass(a)=3;
    
end

Motif3Present=MotifPresent(:,4)+MotifPresent(:,5);
    
%sort filtered sequences for saving
[tempsort,FilteredOrder]=sort(FilteredSequencesCount,'descend');    
    
% Save the data
% Data is saved as tab separated text, with a header 

[filehandle, message]=fopen([OutputFilename '.txt'], 'w');

if(filehandle==-1)
    error(['Error opening file ' OutputFilename ': ' message]);
end    
    

% write header
fprintf(filehandle, 'Sequence\tCount\tFrequency\tLibrary\tATGTT\n');


% write data

for(b=1:1:TotalFilteredNumber)
    
    fprintf(filehandle, '%s\t%u\t%g\t%u\t%u\n', FilteredSequences{FilteredOrder(b)},FilteredSequencesCount(FilteredOrder(b)),FilteredSequencesFreq(FilteredOrder(b)),LibraryClass(FilteredOrder(b)),Motif3Present(FilteredOrder(b)));
        
end

fclose(filehandle);


%% Save the data (Fasta)
% Data is saved as Fasta 

[filehandle, message]=fopen([OutputFilename '.fas'], 'w');

if(filehandle==-1)
    error(['Error opening file ' OutputFilename ': ' message]);
end    
    
% write data

for(b=1:1:TotalFilteredNumber)
    
    fprintf(filehandle, '>No%u_Cnt_%u_Freq_%2.2f_Cl%u\n%s\n', b, FilteredSequencesCount(FilteredOrder(b)),FilteredSequencesFreq(FilteredOrder(b)), LibraryClass(FilteredOrder(b)), FilteredSequences{FilteredOrder(b)});
        
end

fclose(filehandle);

%% Save Matlab file for further analysis

for(b=1:1:TotalFilteredNumber)
    temp=FilteredSequences{FilteredOrder(b)};
    SortedSeq{b}=temp;
    SortedSeqCore{b}=temp(RandomRegionLimits(FilteredOrder(b),1):RandomRegionLimits(FilteredOrder(b),2));
    EnrichmentPercentage(b)=FilteredSequencesFreq(FilteredOrder(b));
    EnrichmentCount(b)=FilteredSequencesCount(FilteredOrder(b));
end

save([OutputFilename '.mat'], 'SortedSeq','SortedSeqCore','EnrichmentPercentage', 'EnrichmentCount');
    
%disp('end')
