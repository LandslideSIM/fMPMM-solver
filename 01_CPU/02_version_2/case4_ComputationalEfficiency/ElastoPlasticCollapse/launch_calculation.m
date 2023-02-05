clear,close

% delete('C:\Users\agnes\Desktop\MPM\work\technical_note\case3_ElastoPlasticCollapse_new\GIMPM_solution\data\*')
% delete('C:\Users\agnes\Desktop\MPM\work\technical_note\case3_ElastoPlasticCollapse_new\GIMPM_solution_iterative\data\*')
% 
% cd 'C:\Users\agnes\Desktop\MPM\work\technical_note\case3_ElastoPlasticCollapse_new\GIMPM_solution_iterative'
% run UNIL_2DGIMP_EP__collapse_iterative.m
% cd 'C:\Users\agnes\Desktop\MPM\work\technical_note\case3_ElastoPlasticCollapse_new\GIMPM_solution'
% run UNIL_2DGIMP_EP__collapse.m

cd 'C:\Users\Manu\Desktop\case3_ElastoPlasticCollapse_new'
run postprocessing.m
run postprocessing_iterative.m