function seqanalyse_sso
%% analysis
% adjusted for Seung Soo, ignores first 5 nucleotides
DEBUG=0;

Filename='d:/Seung Soo/S11.fastq';
OutputFilename='S11_raw';
TopSequencesToSave=1000;

if(~DEBUG)
    tic
    [Header, Sequences, Quality] = fastqread(Filename);
    toc
    NumberOfSequences=length(Sequences);   
else    
    NumberOfSequences=200000;
    tic
    [Header, Sequences, Quality] = fastqread(Filename,'Blockread',[1 NumberOfSequences]);
    toc
end

%remove leading first five bases for Seung Soo in an inefficient loop,
%because cell arrays suck!
for(a=1:1:NumberOfSequences)
    tempstr=Sequences{a};
    tempstr(1:5)=[];
    Sequences(a)={tempstr};
end

SequenceIndex=zeros(NumberOfSequences,1); %Gives an index to each sequence

SequenceCount=zeros(NumberOfSequences,1); %Counts how many times a sequence of with this index exists

index=1;
CurrentPosition=1;

tic
SequencesSorted=sort(Sequences); %Sequence sorting is fast in matlab and greatly accelerates runtime
toc

tic

while(1)
      
    if(CurrentPosition>=NumberOfSequences-1)
        break;
    end
    
    SequenceIndex(CurrentPosition)=index;
    SequenceCount(index)=1;
    
    for(a=CurrentPosition+1:1:NumberOfSequences)
        
        if(SequenceIndex(a)==0)
            
            if(strcmp(SequencesSorted{CurrentPosition},SequencesSorted{a})==1)
                
                SequenceIndex(a)=index;
                SequenceCount(index)=SequenceCount(index)+1;
                
            else % Since Sequences are sorted, stop after first mismatch.
                
                break;
                
            end
            
        end
    
    end
    
    index=index+1;
    CurrentPosition=a;
        
end

toc

frequencies=double(SequenceCount)/double(NumberOfSequences);

[temp,order]=sort(SequenceCount,'descend');

%% Save the data (CSV)
% Data is saved as tab separated text, with a header 

frequencies=frequencies*100; % convert everything into percentages

[filehandle, message]=fopen([OutputFilename '.txt'], 'w');

if(filehandle==-1)
    error(['Error opening file ' OutputFilename ': ' message]);
end    
    

% write header
fprintf(filehandle, 'Sequence\tCount\tFrequency\n');


% write data

for(b=1:1:TopSequencesToSave)
    
    fprintf(filehandle, '%s\t%u\t%g\n', SequencesSorted{order(b)},SequenceCount(order(b)),frequencies(order(b)));
        
end

fclose(filehandle);


%% Save the data (Fasta)
% Data is saved as Fasta 

[filehandle, message]=fopen([OutputFilename '.fas'], 'w');

if(filehandle==-1)
    error(['Error opening file ' OutputFilename ': ' message]);
end    
    
% write data

for(b=1:1:TopSequencesToSave)
    
    fprintf(filehandle, '>SequenceNo_%u_Abundancy_%u_FreqPercentage_%g\n%s\n', b, SequenceCount(order(b)),frequencies(order(b)), SequencesSorted{order(b)});
        
end

fclose(filehandle);


%% save raw data in matlab file
%
save([OutputFilename '.mat'], 'SequencesSorted','SequenceCount','frequencies', 'SequenceIndex','order');

if(DEBUG) %Breakpoint if in Debug mode 
    DEBUG=DEBUG;
end