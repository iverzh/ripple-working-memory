
close all 
clc
clear

%%

subj_list_full = {'P41CS_R1', 'P41CS_R2', 'P42CS_R1', 'P42CS_R2', 'P43CS_R1', 'P43CS_R2', 'P44CS_R1', 'P47CS_R1', 'P47CS_R2', ...
                  'P49CS_R1', 'P49CS_R2', 'P51CS_R1', 'P51CS_R2', 'P53CS_R1', 'P53CS_R2', 'P54CS_R1', 'P54CS_R2', ...
                  'P55CS_R1', 'P55CS_R2', 'P56CS_R1', 'P56CS_R2', 'P57CS_R1', 'P57CS_R2', 'P58CS_R1', 'P60CS_R1', ...
                  'P62CS_R1', 'P62CS_R2'};

matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';

recordingState = 'wake';
location = 'NC';
modifier = '';
tag = [recordingState,'_',location,'_',modifier];

figure;
for subj = 1:length(subj_list_full)
    
    subject = subj_list_full{subj};
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    load(fullfile(matExportFolder, filename))




end























