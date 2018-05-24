
function json_to_list(json_obj) {
    var result_arr = [];
    for (single_entry in json_obj) {
        var sub_arr = Object.values(single_entry.fields)
        result_arr.push(sub_arr);
    }
    return result_arr;
}