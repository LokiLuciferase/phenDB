
function resultList_to_datatable_input(json_obj) {
    var all_data = [];
    var titles = [];

    var all_keys = Object.keys(json_obj[0].fields);
    for(let i=0; i<all_keys.length; i++) {
        let dicTmp = {title : all_keys[i]};
        titles.push(dicTmp);
    }

    for (let i=0; i<json_obj.length; i++) {
        var sub_arr = Object.values(json_obj[i].fields);
        all_data.push(sub_arr);
    }

    return [titles, all_data];
}

//WARNING: This code could break if the Django model is changed
function all_models_to_list(json_obj) {
    var model_names = [];
    var types = [];
    var model_descriptions = [];
    var model_train_dates = [];
    for(let i=0; i<json_obj.length; i++) {
        var row = Object.values(json_obj[i].fields);
        model_names.push(row[0]);
        types.push(row[1]);
        model_descriptions.push(row[2]);
        model_train_dates.push(row[3]);
    }
    return [model_names, types, model_descriptions, model_train_dates];
}