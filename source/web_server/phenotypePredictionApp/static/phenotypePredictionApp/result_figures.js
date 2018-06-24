function ResultFigures(resultsListJSTitles, resultsListJSValues, model_names, model_descriptions, bins)
{

    this.resultsListJSTitles = resultsListJSTitles;
    this.resultsListJSValues = resultsListJSValues;
    this.model_names = model_names;
    this.model_descriptions = model_descriptions;
    this.bins = bins;

    this.initialize_result_figures = function() {
        //DATATABLE 1 (Accuracy & Probability info)
        var titles_table1 = convertTitles(["Bin", "Model", "Prediction", "Model_Confidence", "Balanced_Accuracy"]);
        var dataTable1 = this.__initialize_data_table(this.resultsListJSValues, titles_table1, "#trait_prediction_accuracy_table");
        //DATATABLE 2 (Trait prediction)
        var matrix_for_dt2 = resultslist_to_dt2_matrix(this.resultsListJSValues, this.model_names);
        var dataTable2 = this.__initialize_data_table(matrix_for_dt2[0], matrix_for_dt2[1], "#trait_prediction_table");
        //DATATABLE 3 (Trait Count summary)
        calcTraitCounts(this.resultsListJSValues, 0.5, 0.7, models); //TODO: change values
        var dataTable3;
        //FILTER / AUTOCOMPLETE / INFO
        this.__initialize_pica_models_autocomplete(this.model_names, dataTable1);
        this.__initialize_bins_autocomplete(this.bins, dataTable1);
        this.__initialize_pval_cutoff_spinner(dataTable1);
        this.__initialize_accuracy_cutoff_spinner(dataTable1);
    };

    this.__initialize_data_table = function(resultsListJSValues, resultsListJSTitles, identifier) {
        var dataTable = $(identifier).DataTable({
            "lengthMenu": [[50, 100, -1], [50, 100, "All"]],
            data: resultsListJSValues,
            columns: resultsListJSTitles,
            searching: true,
            dom: '<"table_buttons"B>l<"result_table"t><"table_pagination"p>',
            buttons: [
                {
                    extend: 'csv',
                    text: 'Download table as csv'
                },
                {
                    extend: 'excel',
                    text: 'Download table as Excel (.xlsx)'
                }
            ],
        });
        return dataTable;
    };


    //OLD CODE -> TODO: remove as soon as new layout is accepted
    this.__initialize_pica_models_info = function(model_names, model_descriptions) {
        document.getElementById('dt_results_model_filter_info_text').innerHTML = models_to_infotext(model_names, model_descriptions);
        $('#dt_results_model_filter_info_text').puidialog({
            title: "PICA models",
            resizable: false,
            width: "auto",
            responsive: true
        });
    };

    //

    this.__initialize_pica_models_autocomplete = function(model_names, dataTable) {
        $('#dt_results_model_filter').puiautocomplete({
            completeSource: model_names,
            multiple: true,
        });

        $('#dt_results_model_filter').on('focusin focusout keyup', function () {
            var all_items_htmlcoll = this.parentElement.parentElement.getElementsByTagName('li');
            var all_items = Array.prototype.slice.call(all_items_htmlcoll);
            //makes a collection of the raw text, removes empty entries and builds a regex string out of this collection
            var search_expr = all_items.map(x => x.textContent)
                                        .filter(x => x.length > 0)
                                        .map(x => '^' + x + '$')
                                        .join("|");
            dataTable
                .columns(1)
                .search(search_expr, true, false, true)
                .draw();
        });
    },

    this.__initialize_bins_autocomplete = function(bins, dataTable) {
        $('#dt_results_bin_filter').puiautocomplete({
            completeSource: bins,
            multiple: true,
        });
        $('#dt_results_bin_filter').on('focusin focusout keyup', function () {
            var all_items_htmlcoll = this.parentElement.parentElement.getElementsByTagName('li');
            var all_items = Array.prototype.slice.call(all_items_htmlcoll);
            //makes a collection of the raw text, removes empty entries and builds a regex string out of this collection
            var search_expr = all_items.map(x => x.textContent)
                                        .filter(x => x.length > 0)
                                        .map(x => '^' + x + '$')
                                        .join("|");
            dataTable
                .columns(0)
                .search(search_expr, true, false, true)
                .draw();
        });
    };

    this.__initialize_pval_cutoff_spinner = function(dataTable) {
        //custom filtering function
        $.fn.dataTable.ext.search.push(
            function (settings, data, dataIndex) {
                var min = parseFloat($('#dt_results_pval_filter').val(), 10);
                var value = parseFloat(data[3]) || 0;

                if (isNaN(min) ||
                    ( min <= value)) {
                    return true;
                }
                return false;
            }
        );
        $('#dt_results_pval_filter').keyup(function () {
            dataTable.draw();
        });
    };

    this.__initialize_accuracy_cutoff_spinner = function(dataTable) {
        //custom filtering function
        $.fn.dataTable.ext.search.push(
            function (settings, data, dataIndex) {
                var min = parseFloat($('#dt_results_accuracy_filter').val(), 10);
                var value = parseFloat(data[4]) || 0;

                if (isNaN(min) ||
                    ( min <= value)) {
                    return true;
                }
                return false;
            }
        );
        $('#dt_results_accuracy_filter').keyup(function () {
            dataTable.draw();
        });
    };
}