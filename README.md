HiggsAnalysis-CombinedLimit
===========================

[Manual to run combine](https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit#How_to_run_the_tool)

### Recipe to check out RazorCMS version of combine
```
cmsrel CMSSW_6_1_2
cd CMSSW_6_1_2/src/
cmsenv
git clone https://github.com/RazorCMS/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git pull origin razor1dpdf
scramv1 b clean; scramv1 b

### Recipe for a successful merge with conflicts
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git comb
cd comb/
git checkout master
git merge [branch_to_merge]
git checkout -b [branch_to_merge]
git pull git@github.com:[user_to_merge]/HiggsAnalysis-CombinedLimit.git [branch_to_merge]
[resolve conflicts]
git add [files where conflicts were resolved]
git commit #(should contain a header telling something about merging [branch_to_merge])
git push origin master
```
