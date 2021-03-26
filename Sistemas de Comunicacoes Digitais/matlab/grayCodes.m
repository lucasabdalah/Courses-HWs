function Seq=grayCodes(nbit)
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
    base = 2;
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
    