        //var resultsList_arr = resultList_to_datatable_input({{ resultsList|jsonify | safe}});
        //var resultsListJSTitles = resultsList_arr[0];
        //var resultsListJSValues = resultsList_arr[1];

        //var models_arr = all_models_to_list({{ all_models|jsonify | safe }});
        //var model_names = models_arr[0];
        //var model_descriptions = models_arr[1];

        //var dataTable;

        function initialize_result_figures(resultsListJSTitles, resultsListJSValues, model_names, model_descriptions, bins) {
            var dataTable =  __initialize_data_table(resultsListJSValues, resultsListJSTitles);
            __initialize_pica_models_info(model_names, model_descriptions);
            __initialize_pica_models_autocomplete(model_names, dataTable);
            __initialize_bins_autocomplete(bins, dataTable);
        }

        function __initialize_data_table(resultsListJSValues, resultsListJSTitles) {
            var dataTable = $('#dt_results_table').DataTable( {
                            data: resultsListJSValues,
                            columns: resultsListJSTitles,
                            searching: true
            } );
            return dataTable;
        }

        function __initialize_pica_models_info(model_names, model_descriptions) {
            document.getElementById('dt_results_model_filter_info_text').innerHTML = models_to_infotext(model_names, model_descriptions);
            $('#dt_results_model_filter_info_text').puidialog({
                title: "PICA models",
                resizable: false,
                width: "auto",
                responsive: true
            });
        }

        function __initialize_pica_models_autocomplete(model_names, dataTable) {
             $('#dt_results_model_filter').puiautocomplete({
                completeSource: model_names,
                multiple: true,
            });
            $('#dt_results_model_filter_info').puibutton({
                icon: 'fa-external-link-square',
                click: function() {
                    $('#dt_results_model_filter_info_text').puidialog('show');
                }
            });

            $('#dt_results_model_filter').on('focusin focusout keyup', function() {
                var all_items_htmlcoll = this.parentElement.parentElement.getElementsByTagName('li');
                var all_items = Array.prototype.slice.call(all_items_htmlcoll);
                var search_expr = all_items.map(x => x.innerText).filter(x => x.length > 0).map(x => '^' + x + '$').join("|");
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
                var search_expr = all_items.map(x => x.innerText).filter(x => x.length > 0).map(x => '^' + x + '$').join("|");
                dataTable
                    .columns(0)
                    .search(search_expr, true, false, true)
                    .draw();
            });
        }

        function __initialize_verdict_cutoff_spinner(dataTable) {
            $.fn.dataTable.ext.search.push(
                function( settings, data, dataIndex ) {
                    var min = parseInt( $('#min').val(), 10 );
                    var verdict = parseFloat( data[3] ) || 0;

                    if  ( isNaN( min ) ||
                        ( isNaN( min ) && verdict <= max ) ||
                        ( min <= age) )
                    {
                        return true;
                    }
                    return false;
                }
            );
            $('#dt_results_verdict_filter').keyup(function() {
                dataTable.draw();
            })
        }
