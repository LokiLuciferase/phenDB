var all_keys;

function json_to_datatable_input(json_obj) {
    var all_data = [];
    var titles = []

    all_keys = Object.keys(json_obj[0].fields);
    for(i=0; i<all_keys.length; i++) {
        console.log(all_keys[i])
        let dicTmp = {title : all_keys[i]};
        titles.push(dicTmp);
    }

    for (i=0; i<json_obj.length; i++) {
        var sub_arr = Object.values(json_obj[i].fields);
        all_data.push(sub_arr);
    }

    return [titles, all_data];
}