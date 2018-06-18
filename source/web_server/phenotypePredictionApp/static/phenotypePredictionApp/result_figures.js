    function initialize_result_figures(resultsListJSTitles, resultsListJSValues, model_names, model_descriptions, bins) {
        var titles_table1 = __convertTitles(["Bin", "Model", "Prediction", "Model_Confidence", "Balanced_Accuracy"]);
        var dataTable1 =  __initialize_data_table(resultsListJSValues, titles_table1, "#trait_prediction_accuracy_table");
        var matrix_for_dt2 = resultslist_to_dt2_matrix(resultsListJSValues, model_names);
        var dataTable2 = __initialize_data_table(matrix_for_dt2[0], matrix_for_dt2[1], "#trait_prediction_table"); //unshift adds empty item to the beginning of the array
        __initialize_pica_models_autocomplete(model_names, dataTable1);
        __initialize_bins_autocomplete(bins, dataTable1);
        __initialize_pval_cutoff_spinner(dataTable1);
        __initialize_accuracy_cutoff_spinner(dataTable1);
    }

    function __initialize_data_table(resultsListJSValues, resultsListJSTitles, identifier) {
        window.resultsListJSValues = resultsListJSValues;
        window.resultsListJSTitles = resultsListJSTitles;
        var dataTable = $(identifier).DataTable( {
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
        } );
        return dataTable;
    }


    //OLD CODE -> TODO: remove as soon as new layout is accepted
    function __initialize_pica_models_info(model_names, model_descriptions) {
        document.getElementById('dt_results_model_filter_info_text').innerHTML = models_to_infotext(model_names, model_descriptions);
        $('#dt_results_model_filter_info_text').puidialog({
            title: "PICA models",
            resizable: false,
            width: "auto",
            responsive: true
        });
    }
    //

    function __initialize_pica_models_autocomplete(model_names, dataTable) {
         $('#dt_results_model_filter').puiautocomplete({
            completeSource: model_names,
            multiple: true,
        });

        $('#dt_results_model_filter').on('focusin focusout keyup', function() {
            var all_items_htmlcoll = this.parentElement.parentElement.getElementsByTagName('li');
            var all_items = Array.prototype.slice.call(all_items_htmlcoll);
            //makes a collection of the raw text, removes empty entries and builds a regex string out of this collection
            var search_expr = all_items.map(x => x.textContent).filter(x => x.length > 0).map(x => '^' + x + '$').join("|");
            dataTable
                .columns(1)
                .search(search_expr, true, false, true)
                .draw();
        });
    }

    function __initialize_bins_autocomplete(bins, dataTable) {
         $('#dt_results_bin_filter').puiautocomplete({
            completeSource: bins,
            multiple: true,
        });
        $('#dt_results_bin_filter').on('focusin focusout keyup', function() {
            var all_items_htmlcoll = this.parentElement.parentElement.getElementsByTagName('li');
            var all_items = Array.prototype.slice.call(all_items_htmlcoll);
            //makes a collection of the raw text, removes empty entries and builds a regex string out of this collection
            var search_expr = all_items.map(x => x.textContent).filter(x => x.length > 0).map(x => '^' + x + '$').join("|");
            dataTable
                .columns(0)
                .search(search_expr, true, false, true)
                .draw();
        });
    }

    function __initialize_pval_cutoff_spinner(dataTable) {
        //custom filtering function
        $.fn.dataTable.ext.search.push(
            function( settings, data, dataIndex ) {
                var min = parseFloat( $('#dt_results_pval_filter').val(), 10 );
                var value = parseFloat( data[3] ) || 0;

                if  ( isNaN( min ) ||
                    ( min <= value) )
                {
                    return true;
                }
                return false;
            }
        );
        $('#dt_results_pval_filter').keyup(function() {
            dataTable.draw();
        });
    }

    function __initialize_accuracy_cutoff_spinner(dataTable) {
        //custom filtering function
        $.fn.dataTable.ext.search.push(
            function( settings, data, dataIndex ) {
                var min = parseFloat( $('#dt_results_accuracy_filter').val(), 10 );
                var value = parseFloat( data[4] ) || 0;

                if  ( isNaN( min ) ||
                    ( min <= value) )
                {
                    return true;
                }
                return false;
            }
        );
        $('#dt_results_accuracy_filter').keyup(function() {
            dataTable.draw();
        });
    }
