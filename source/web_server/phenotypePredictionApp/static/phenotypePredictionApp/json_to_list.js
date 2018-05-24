
function json_to_list(json_obj) {
    var result_arr = [];
    for (single_entry in jsob_obj) {
        var sub_arr = [];
        for (fieldKey in single_entry.fields) {
            sub_arr.push(single_entry.fields[fieldKey]);
        }
        result_arr.push(sub_arr);
    }
    return result_arr;
}