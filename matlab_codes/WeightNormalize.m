function P=WeightNormalize(P,R)
%function to normalize prediction score from REMAP
%Normalization weight 5.25 by default
%weight derived based on ZINC_ChEMBL dataset with negative examples
%P: prediction score matrix
%R: observed pair adjacent matrix
%bins=[0 to 1.10], bin width=0.05

wt=5.25; %weight on negative pair counts
bin_centers=0.025:0.05:1.1;
bin_weights=zeros(size(bin_centers));
for bin_idx=1:length(bin_centers)
    bin_positive=sum(R((P>(bin_centers(bin_idx)-0.025)) & P<((bin_centers(bin_idx)+0.025))));
    bin_total=length(R((P>(bin_centers(bin_idx)-0.025)) & P<((bin_centers(bin_idx)+0.025))));
    bin_negative=bin_total-bin_positive;
    bin_weights(1,bin_idx)=(bin_positive/(bin_positive+5.25*bin_negative)); %adjusted score for the bin
end

[n,m]=size(P);
total=n*m; %total number of pairs for progress report
count=0;
tenpercent=round(total/10);
for i=1:n
    for j=1:m
        score=P(i,j);
        for bin_idx=1:length(bin_centers)
            if (score>(bin_centers(bin_idx)-0.025) && score<(bin_centers(bin_idx)+0.025) )
                P(i,j)=bin_weights(bin_idx); %new score (adjusted)
                count=count+1;
                if mod(count,tenpercent)==0
                     disp(['Normalization... ' num2str(10*count/tenpercent) ' percent done']);
                end
                break
            end
        end
    end
end

end