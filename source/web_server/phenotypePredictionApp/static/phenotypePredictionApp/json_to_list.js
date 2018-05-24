test;

function json_to_list(json_obj) {
    var result_arr = [];
    for (i=0; i<json_obj.size(); i++) {
        test = single_entry[i].fields
        var sub_arr = Object.values(single_entry[i].fields);
        result_arr.push(sub_arr);
    }
    return result_arr;
}