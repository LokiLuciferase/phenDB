
function resultList_to_datatable_input(json_obj) {
    var all_data = [];
    var titles = [];

    var all_keys = Object.keys(json_obj[0].fields);
    titles = convertTitles(all_keys);

    for (var i=0; i<json_obj.length; i++) {
        var sub_arr = Object.values(json_obj[i].fields);
        all_data.push(sub_arr);
    }
    __roundNumbers(all_data, [3,4], 2);
    return [titles, all_data];
}

function convertTitles(titles_arr) {
    var titles_new = [];
    for(var i=0; i<titles_arr.length; i++) {
        var dicTmp = {title : titles_arr[i]};
        titles_new.push(dicTmp);
    }
    return titles_new;
}

function __roundNumbers(all_data, columns, precision) {
    for(var i=0; i<all_data.length; i++) {
        for(var a=0; a<columns.length; a++) {
            var col = columns[a];
            var precision_factor = 10*precision;
            var value_rounded = Math.round(all_data[i][col] * precision_factor)/precision_factor;
            var value__rounded_fixed_decimal = Number.parseFloat(value_rounded).toFixed(precision);
            all_data[i][col] = value__rounded_fixed_decimal;
        }
    }
}

//WARNING: This code could break if the Django model is changed
//TODO: model description may bot be necessary anymore
function all_models_to_list(json_obj) {
    var model_names = [];
    var model_descriptions = [];
    for(var i=0; i<json_obj.length; i++) {
        var row = Object.values(json_obj[i].fields);
        if(model_names.indexOf(row[0]) != -1) continue;
        model_names.push(row[0]);
        model_descriptions.push(row[2]);
    }
    return [model_names, model_descriptions];
}

//OLD CODE -> TODO: remove as soon as new layout is accepted
function models_to_infotext(titles, descriptions) {
    var infotable = "<table>";
    for(var i=0; i<titles.length; i+=2) {
        infotable += "<tr>";
        for(var a=0; a<2; a++) {
            infotable += "<td class='table_info_box'>";
            infotable += titles[i+a];
            infotable += "</td>";
            infotable += "<td class='table_info_box'>";
            if(descriptions[i+a].length == 0) infotable += "no description available";
            else infotable += "description: " + descriptions[i+a];
            infotable += "</td>";
        }
        infotable += "</tr>";
    }
    infotable += "</table>";
    return infotable;
}
//

function getAllBins(resultsListJSValues) {
    var bins = resultsListJSValues.map(x => x[0]);
    var unique_bins = bins.filter(function(value, position) {
        return bins.indexOf(value) == position;
    });
    return unique_bins;
}

/*
function resultslist_to_dt2_matrix(resultsListJSValues) {
    var bin_column = window.database_structure.bin_table.bin_column;
    var model_column = window.database_structure.bin_table.model_column;

    function unique(value, index, self) {
        return self.indexOf(value) === index;
    }

    var bins = resultsListJSValues
                        .map(x => x[bin_column])
                        .filter(unique);

    var models = resultsListJSValues
                        .map(x => x[model_column])
                        .filter(unique);
} */


function resultslist_to_dt2_matrix(resultsListJSValues, models) {
    resultsListJSValues.sort(function (a, b) {
        if (a[0] < b[0]) return -1;
        else if (a[0] > b[0]) return 1;
        else return 0;
    });

    var dt2_matrix = [];
    var bin_prev = resultsListJSValues[0][0];
    var row_index = 0;
    var row = [];
    row[0] = bin_prev;
    for(var i=0; i<resultsListJSValues.length; i++) {
        var bin_tmp = resultsListJSValues[i][0];
        var verdict = resultsListJSValues[i][2];
        var model = resultsListJSValues[i][1];
        var model_index = models.indexOf(model);
        if(model_index == -1) throw "model not in model set";
        row[model_index+1] = verdict;
        if(bin_prev != bin_tmp || i == resultsListJSValues.length-1) {
            row = __initializeEmptySlots(row);
            dt2_matrix[row_index] = row;
            row = [];
            row[0] = bin_tmp;
            ++row_index;
        }
        bin_prev = bin_tmp;
    }

    var model_cp = models.slice();
    model_cp.unshift("bins");
    var titles = convertTitles(model_cp);
    return [dt2_matrix, titles];
}

function __initializeEmptySlots(arr) {
    for(var i=0; i<arr.length; i++) {
        if (arr[i] == undefined) arr[i] = "NA";
    }
    return arr;
}


function calcTraitCounts(resultsListJSValues, confidenceCutoff, accuracyCutoff, models) {
    var false_dir = {};
    var true_dir = {};
    var nd_dir = {};
    var nc_dir = {};

    var model_column = window.database_structure.bin_table.model_column;
    var prediction_column = window.database_structure.bin_table.prediction_column;
    var confidence_column = window.database_structure.bin_table.confidence_column;
    var accuracy_column = window.database_structure.bin_table.accuracy_column;

    for(var i=0; i<resultsListJSValues.length; i++) {
        var model = resultsListJSValues[i][model_column];
        var prediction = resultsListJSValues[i][prediction_column];
        var confidence = resultsListJSValues[i][confidence_column];
        var accuracy = resultsListJSValues[i][accuracy_column];
        if(confidence < confidenceCutoff || accuracy < accuracyCutoff) nd_dir[model] = (nd_dir[model] || 0) +1;
        else if(prediction === false) false_dir[model] = (false_dir[model] || 0) + 1;
        else if(prediction === true) true_dir[model] = (true_dir[model] || 0) +1;
    }

    var result_matrix = [];

    for(var i=0; i<models.length; i++) {
        var row = [];
        var model_tmp = models[i];

        row.push(model_tmp);

        if (true_dir[model_tmp] == undefined) row.push(0);
        else row.push(true_dir[model_tmp]);

        if(false_dir[model_tmp] == undefined) row.push(0);
        else row.push(false_dir[model_tmp]);

        if(nd_dir[model_tmp] == undefined) row.push(0);
        else row.push(nd_dir[model_tmp]);

        result_matrix.push(row);
    }

    var titles = convertTitles(["", "YES", "NO", "ND"]);
    return [result_matrix, titles];
}