function symb= mapping_MQAM(M,bits)


    switch M
    
    % Se estivermos trabalhando com 4-QAM
    case 4 
        switch bits
            case [0000]
                symb = 
            
            
        end
    
    % Se estivermos trabalhando com 16-QAM
    case 16
    
    % Se estivermos trabalhando com 64-QAM
    case 64
    
    end

    % % M-QAM mapping
    % if bits==[0]
    %     symb= -1;
    % end
    % if bits==[1]
    %     symb= 1;
    % end
    
% To get Gray-code back to binary it needs only to sort()



https://www.mathworks.com/matlabcentral/fileexchange/62403-gray-code-for-all-base-numbers
function Seq=grayCodes(base,nbit)
    %% Gray Code Generation
    % 
    % # Er.Abbas Manthiri S
    % # abbasmanthiribe@gmail.com
    % # Matlab 2014a
    % # 04-04-2017
    % 
    % Inputs
    % nbit=single value;
    % base=single value;
    % Output
    % Seqence
    SeqLength=base^nbit;
    Seq=zeros(SeqLength,nbit);
    for i=1:nbit
        RepetaionRange=base^(nbit-i+1)/base;
        IncrementSeq=0;
        sequenceElement=1;
        while(IncrementSeq<SeqLength)
            end_point=IncrementSeq+RepetaionRange;
            start_point=IncrementSeq+1;
            Seq(start_point:end_point,i)=sequenceElement-1;
            IncrementSeq=IncrementSeq+RepetaionRange;
            sequenceElement=sequenceElement+1;
            if sequenceElement>base
                sequenceElement=1;
            end
        end
    end
    