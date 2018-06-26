window.database_structure = {
    bin_table : {
        bin : {db_name : "bin", ui_name: "Bin"},
        model : {db_name : "model", ui_name: "Model"},
        prediction : {db_name: "verdict", ui_name:"Prediction"},
        confidence : {db_name: "pica_pval", ui_name:"Confidence"},
        accuracy : {db_name: "accuracy", ui_name:"Accuracy"},
        nc_masked: {db_name: "nc_masked", ui_name:""}, //won't be directly displayed in UI, needed for logic
    },
};