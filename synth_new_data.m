projectRoot  = "C:\Users\yariv\OneDrive - Afeka College Of Engineering\BSc_EEE\Final_BSc_Project\Part B\SM_Correlation";
customersDir = fullfile(projectRoot, "data", "virtual_customers");
outXlsx      = fullfile(projectRoot, "data", "transformer_profile.xlsx");

dir(fullfile(customersDir,"*.xlsx"))   % בדיקה שחייב להחזיר קבצים
generate_customers(120, 100);
generate_transformer_profile(customersDir, outXlsx);
