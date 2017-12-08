function merged_table = make_MERRA2_table_TCO_PWV_column_multi_years()
%years = 2006:1:2015;
%years = 2016;
years = 2017;
N = size(years);
merged_table = table; 
for j = 1:1:N(2)
   %MERRA2_column = make_MERRA2_table_TCO_PWV_column(years(j));
   MERRA2_column = make_MERRA2_table_TCO_PWV_partial_column(years(j));
   merged_table = [merged_table;MERRA2_column];
end