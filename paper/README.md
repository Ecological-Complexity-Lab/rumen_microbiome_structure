# Recommendation by the Subject Editor:   
  
We have now received the feedback from two reviewers and they are positive about the study, and I agree with them. I enjoyed reading your manuscript. I found it interesting and well-written, and the analyses robust. Yet, the reviewers provided a number of relevant suggestions that need to be addressed before considering the manuscript for publication. I believe that addressing such comments and suggestions will help to further improve the quality and scope of your manuscript. I hope you find them useful.   

- Thank you for that positive feedback. We have addressed the reviewers' comments, which indeed improved the quality of the manuscript. We present the point-by-point answers below.


# Reviewer: 1   
  
Comments to the Author   
This paper presents an interesting, multilayer co-occurrence network analysis for rumen microbiomes inhabiting cows in different farms. The analysis is technically and statistically correct, I mostly ave concerns about the principles and the concepts. But overall it is a nice and useful study.   
  
lines 37-38: these three components of "community membership" are much older in the ecological literature, dating back to the mid-seventies, please check, for example, this classical piece of work:   
https://www.degruyter.com/document/doi/10.1515/9781400857081.297/pdf?licenseType=restricted   
- [ ] We have  aded this paper to the citations


line 43: You may have co-occurrence without any interaction. For example, f there is a strong abiotic factor in the environment (e.g. iron, salt), organisms preferring it may co-exist without any interaction. The signatue of this in network terms is a large clique (if we speak about positive associations).   
For an example for higher organisms, an excellent study:   
https://link.springer.com/article/10.1556/ComEc.11.2010.1.14   

- [ ] Yes, this is true. We have rewritten this paragraph and cited this paper. #TODO  

line 44-45: "local" and "proximity" should be defined better, especially if we speak about microbes. If there is metabolic coupling between two microbes, what is proximity? How long is it?   

- [ ] Define #TODO 

line 47: I am not sure that the spatial scale of physical filtering is neccessarily much larger than the scale of competition. Can you justify this?   

- [ ] We have rewritten this paragraph, considering this comment #TODO 
  
line 57: complex relative to what?  It is too easy to say that everything is just complex.   
- [ ] We removed these words #TODO

line 67: the "core" may not only be a function of dominance or frequency, see this paper for quantifying "surprising hubs" in similar networks:   
https://www.nature.com/articles/srep15920   

- [ ] This is true, but we cannot consider this in our system. We included this idea in the discussion #TODO 

In general, the discussion on interaction types (i.e. metabolic, symbiotic etc...) and especially the signs (e.g. negative) and sign combinations (e.g. minus-minus) should be more elaborated. These are crucial for undertanding how are co-occurences and interactions really related to each other.   

- [ ] Add a few sentences to the introduction #TODO 
  
line 72-73: Co-occurrence is always direct, I guess. Statistically it is detected between nodes i and j. Whether the real interaction between them, resulting in the statistical co-occurrence, is direct or indirect is another issue. Please clarify this.   
  - [ ] By 'indirect' we regard the processes that generate cooccurrence or those that emanate from them. We have reworded this sentence #TODO 

H1 makes only sense with positive co-occurrences. Please clarify this. Negative associations are not transitive, instead, they lead to positive associations in two steps (i - j - k).   
- [ ] Yes, this is true we have calrified this #TODO 
  
like 119-120: I am sure about the cows, but is it for sure abut the microbes, too?   
  - [ ] Yes, this is true also for the microbes in this system because they are bound to the cow rumen.
  
line 167: again, please make it clear that you speak about positive co-occurrences.   
- [ ] We have clarified this here too #TODO 

Figure 2 may be related to the total number of microbe individuals per species (abundance)? The real signal would be a wide-spread yet not very abundant microbe.   
- [ ] The figure presents the occurrence (disregarding abundances). The claim may be true, but we are not able to test it, even if we used reads (reads are not a good measure of absolute abundance). We have added the following sentence to the legend: In this study occurrence does not consider microbe abundance. #TODO 
  
line 331-333: this is unclear to me, can you please clarify? This may depend on the similarity of the regional species pool that, I assume, is not very different.  

- Yes, there is a dependency on the regional pool (i.e., union of microbes across farms). The fact that there are shared microbes between farms is beneficial for this hypothesis.  Here is an elaborated explanation of the hypothesis: Lets take an extreme scenario in which the processes that determine co-occurrence are the same in all farms. In that case, microbes are expected to co-occur with the same partners in all farms. The interlayer links will connect microbes across all farms effectively creating a network that is highly connected within and between farms. In the other extreme scenario, the processes that determine cooccurrence are very different between farms. In that case, microbes will co-occur with the same partners in different farms to a very low extent. In that scenario the multilayer network will be highly connected within farms and hardly at all between farms. The modularity (i.e., clustering) algorithm will pick up such scenarios. In the first scenario there will be no obvious separation of farms to clusters because random walks will be performed equaliy between and within farms. In the latter scenario, random walks will tend to concentrate on the dense areas (in the farms), and hardly move between farms, creating a clustered network in which farms are separated to clusters.

The multilayer approach sounds very good but it would be more logical to me if applied to classical interactions (local rules causing modules, yes). For co-occurrence, as it was argued, proximity is important, generating a strong modularity signal, right? 
- [ ]  Following the explanation above, the multilayer approach is necessary to test this hypothesis. Generating a network based solely on cooccurrence across farms will wipe out the farms from the network impeding testing this hypothesis. We emphasize this point #TODO
   
  
It is a pity that the two clusters of cows (Holstein vs red) also receive different diet. I guess many differences are more because of the diet and not because of genetics, but we cnnot see this.   
  - [ ] We agree, it is a disadvantage. However, the diets are equal in their nutrients. <span style="background:#fff88f">Consult Itzik</span>
  
line 397: (Saravia et al. 2022)   
- [ ] We fixed this #TODO
  
line 404-405: yes, stool analysis could be different but also different segments of the digestive tract. There is a huge variability along the oral cavity - stool axis.   
 - [ ] True, we added this comment to the sentence #TODO 
  
Back to negative associations. I think thee are technically harder to handle but might be even more informative. Positive associations can rise easily (e.g. a microelement patch in the environment), while negative ones are really informative and, probably infer real interactions more often. Line 142 mentions negative associations but this remains undiscussed.   
- [ ] We agree, yet these are not the focus of this study and require other hypotheses than those we have posed. Nevertheless they do require at least some discussion as suggested, and we added this to the discussion. #TODO 
  
The "gas emission mitiation" line of thinking appears first in the abstract and second (and last time) at the very end. I think it can be omitted: this idea is not really elaborated and also not really needed for this paper. It is a possible perspective but there are many other potential applications. If it remains there, the Reader has million of questions to answer. For example, a much better analysis of the diet - microbes - co-occurrence - system dynamics - gas emission chain of factors.   
  - [ ] Yes, we agree. Wer have removed this #TODO 

As a general comment, the Authors adopt too easily the ecological concepts from higher organisms to microbes, I think. Ineractions (if thee are real interactions) may be mch more like ad hoc and transitory, the orgnisms have flexible genetic and metabolic systems and they are driven much more by abiotic graidents/signals - in case of microorganisms. I think this could be discussed, even in a table with similarities and dissimilarities between eukaryotes and prokaryotes. Maybe in a next paper.   
- [ ] Yes, we agree. We now mention this in the introduction and discussion #TODO 


# Reviewer: 2   
  
Comments to the Author   
Understanding the impact of microbial community structure on ecosystem function is crucial in modern biology. However, microbial communities are frequently characterized by high diversity, making it challenging to discern their interactions. In this study, the authors tackle this challenge by investigating the signatures of microbial co-occurrence networks across diverse spatial scales. Their analysis encompasses a range of network properties, augmented by simulation and permutation tests.   
I found the study to be engaging and well-written with thoughtful analyses. I recommend accepting this manuscript for publication in Ecography pending the incorporation of my feedback outlined below.   

- Thank you for this positive feedback. We have addressed the comments as described below, which significantly improved the quality of the work.

## Major comments   
I would argue a uni-modular pattern observed in samples from 7 farms in four European countries is hardly global. I agree the spatial scale is more than regional. Maybe continental? 

- [ ] Indeed, this is true. There is a problem with the terminology. By global we did not mean earth-scale, but simply the opposite of non-local. This terminology is confusing, and we therefore changed it to non-local. #TODO 

What processes may lead to high transitivity may depend on the underlying biology, such as the functional guilds, in the case of microbes, dictated by metabolism and cross-feeding. I am curious about the node identity. One alternative to experimental characterization would be perform the same analysis at the family or phylum level instead of ASVs. I wonder whether the detected transitivity at the ASVs level may result, either from the cross-feeding interdependency between taxa belonging to different functional families, or from the neutral competition between taxa belonging to the same functional families. 

- [ ] We performed this analysis and... #TODO #Geut

What is the detection threshold for the presence/absence of a microbe in cows? Does including only abundant microbes (>0.1%, >1% or > 5%) change the scale at which network signatures is detected?   

- [ ] We need to perform this analysis #TODO #Geut 

Another network properties besides clustering coefficient is network motifs, where fully linked motifs directly imply the trio co-occurrence [Milo et al 2002 10.1126/science.298.5594.824]. It's been shown the trophic networks exhibit matching network pattern between local and regional scales in terms of motifs [Baiser et al 2015 10.1111/oik.02532]. I am not aware of motif analysis in co-occurrence networks and this dataset seems to be well suited.   

- Network motifs are particularly meaningful in directed networks because there are many combinations. Both of these studies were performed on directed networks. Our network is undirected. In that case, motifs are not so useful and resemble clustering coefficient. For instance, the apparent competition and exploitative competition motifs in [Baiser et al 2015 10.1111/oik.02532] will be the same one in an undirected network and resemble an unclosed triangle. Similarly, the omnivory motif will be a triangle regardless of the direction of the arrow. The idea of motifs has occurred to us in an early stage of this research, as in the lab we are very aware of research such as the references mentioned and others (e.g., Stouffer et al 2012  doi:10.1126/science.1216556  and [Saravia et al 2022]) but for the reasons described we did not proceed with analyzing motifs.


## Minor comments   
Line 74. Why would the authors expect non-transitive signature in a co-occurrence network? If stochastic dynamics governs communtiy assembly, for instance, in the case of catipillar gut microbiome [Hammer et al 2017. 10.1073/pnas.1707186114], I would expect no transitivity in the co-occurrence netwroks of such system. Could the authors provide the rationale of this hypothesis?   

- [ ]  We agree that under stochastic processes transitivity is not expected. This statement does not contradict the hypothesis. In light of this and other comments we rewrote this paragraph to better explain our hypothesis.

Line 77-81. A rock-paper-scissor dynamics in a competitive network (a directed graph) is often considered as nontransitivity. This is different the definition of transitivity in co-occurrence graphs. Maybe use another example?   

- [ ] Yes, the terminology is similar. We used another example when rewriting the hypothesis.

Line 92-93 on H3. What processes do the author refer to? Ecological processes? Community processes?   
- [ ] We added examples for processes.

How similar are the rumen communities within and between farms? The max Jaccard similarity  between two farms is ~0.6 according to FigS5. I wonder if this is driven by high variation between farm or high similarity between farms. It would be informative to see the JDs between cows within a farm   

- [ ] We have calculated the Jaccard index within a farm and present the results in... #TODO #Geut 


Line 141-143. If I understand correctly, positive and negative cooccurrence links are represented the same (value 1). Is it the case? If so, are the link properties similar? Does accounting for only positive links alter the findings?   

- [ ] We only used positive interactions. We now specify so #TODO 

Line 148-149. Not sure I understand this rationale. How is this probability calculated? Is it the ASV prevalence in a farm (number of cows with this microbe divided by the number of cows in the same farm)?  

- [ ] We now present the calculation of the probability #TODO #Geut 

Line 159-166. Is the randomized network sampling the microbes with or without replacement?   

- [ ] Without replacement. That is, a microbe can be selected and placed in a cow only once. We now specify that #TODO 

Line 173-175. Can the provide the degree distribution for each farm network and describe the fraction of ASVs kept for CC analysis?   

- [ ] Yes, we provide it now #TODO #Geut 

Line 295. Does "Bacteroidetes (52.7%)" mean 52.7% of 946 microbes?  

- [ ] Yes, it means  52.7% of 946 microbes #Geut 

Line 301. The between-family Jaccard distance is an average metric. How different are the communities between cows within a farm versus between cows from different farm?   

- [ ] We performed this analysis in response to minor comment number 4, above.


Fig3. What's the sample size for each pie? What is the null expectation that the significance is comapred against? The pie chart 

- [ ] We added the sample size for each pie in the figure. The null expectation is explained in the Methods section "Calculating significance using z-scores"


Fig3B. Are the farms ordered according to something?   

- [ ] #Geut 

Fig4D. What are the empty cells? 

- [ ] Empty cells indicate that no module contained a given amount of farms. We added this to the figure caption. #TODO 

FigS2. What is NA? What's the number of ASVs in each farm?  

- [ ] #Geut 

FigS3A. It seems there are two peaks (~150 and ~350). Are these cows from distinct farms or spread out evenly across farms?   

- [ ] #Geut 

FigS4. Sample size?   

- [ ] #Geut 

FigS5. I am not sure how this is calculated. Are the data from different cows from the same farm stacked? Or is the the mean distance between cows from two different farms?   

- In this figure the data are aggregated across cows within each farm. However, in response to minor comment 4, we added an analysis within a farm (see that response)

Fig S6. Which network properties are included in the PCA? CC, groups? How much variation is explained by PC1 and PC2?

- [ ] We added this information to the figure legend. #Geut 
