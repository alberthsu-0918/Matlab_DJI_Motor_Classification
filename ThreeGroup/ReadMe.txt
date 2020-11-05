run 
1. pca_Dajan_3group.m to create 3 type classification model
2. group3_transfer.m to classify wave files under test to type 1(pass) 2(fail) or 3(minor)
3. using pca_Dajan_subModel.m to create minor fail / pass model.
4. Using subGroup_transfer.m to test those wav files classified as 3 in step 2. 
1 is pass, 2 is fail.