classdef (Abstract) TableAndTimeTableUtilities
% Static methods (utility functions) to validate and manipulate time series 
% data stored in tables/timetables
%
% Description:
%
%   This abstract class is a repository for static methods related to the 
%   validation and manipulation of tables/timetables that store time series 
%   data, and the conversion of tables/timetables to conventional 2-D matrices 
%   (univariate data) and 3-D arrays (multivariate data), and vice versa.
%
% Naming Convention:
%
%   Generally, if a given utility method includes the word "Table" in its
%   name, then it is applicable to both tables and timetables. If the name
%   includes the word "TimeTable", then it is applicable to timetables only.
%
% Utility Categories:
%
%   Broadly speaking, there are 3 categories of utilities:
%
%    (1) Consistency and error checking utilities that answer questions: 
%
%        The question is indicated by the name of the utility and begins with
%        the verb "are" or "is" and the answer is optionally returned as a 
%        logical scalar:
%
%        o areTableVariableNamesConsistent
%        o areTimeTableDateTimesConsistent
%        o isTableDataNumericWithOutMissingData
%        o isTableDataSinglePath
%        o isTableFormatValid
%        o isTimeTableRegular
%
%        Utilities that begin with "are" test consistency between two time
%        series, while those that begin with "is" test internal consistency
%        of a single time series.
%
%        When called with no outputs these utilities run silently or throw 
%        errors.
%
%    (2) General validation utilities:
%
%        These utilities perform validation actions and return no outputs such
%        that they run silently or throw errors:
%
%        o ensureTimeSeriesTypeConsistency
%        o validateTablesAndTimeTables
%            Convenience wrapper around single series consistency utilities:
%              o isTableDataNumericWithOutMissingData
%              o isTableFormatValid
%              o isTimeTableRegular
%
%    (3) Time series conversion and metadata update utilities:
%
%        These utilities accept time series of one type and format and convert
%        them to another type and format (e.g., convert matrices and arrays 
%        to tables and timetables, or vice versa), or update metadata of 
%        existing tables and timetables:
%
%        o perPathMatrix2PerVariableTable
%        o perPathArray2PerVariableTable
%        o perVariableTable2PerPathMatrix
%        o perVariableTable2PerPathArray
%        o copyMetaData2Table
%
% Notes:
%
%   o When creating regular timetables, it's better to allow the constructors
%     "timetable" or "array2timetable" to compute the dates rather than 
%     specifying the dates explicitly. In other words, use the N-V pairs 
%     'StartTime' and 'TimeStep' rather than 'RowTimes'. The former syntax 
%     may be thought of as an "implicit" date specification while the latter
%     "explicit". 
%
%     The reason for this recommendation is that when indexing into a timetable 
%     TT to access a subset of the original data, TT.Properties.TimeStep is
%     generally preserved if TT was created with 'StartTime' and 'TimeStep'.
%     In contrast, if TT was created with 'RowTimes' then TT.Properties.TimeStep
%     is only preserved if the subset of TT has two or more rows, otherwise
%     TT.Properties.TimeStep = NaN. 
%
%     A non-NaN TimeStep facilitates the determination of regularity and
%     allows dates into the future and past to be computed from TT.
%

% Copyright 2019 The MathWorks, Inc.

methods (Static)

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = areTableVariableNamesConsistent(TT1,TT2,names)
% True if variables names of tables/timetables are the same
%
% Syntax:
%   I =  areTableVariableNamesConsistent(TT1,TT2,Names)
%   areTableVariableNamesConsistent(TT1,TT2,Names)
%
% Description:
%   Given tables/timetables TT1 and TT2, determine if the variable names of
%   TT1 and TT2 are the same (i.e., TT1.Properties.VariableNames is the same
%   as TT2.Properties.VariableNames).
%
%   When called with an output, true is returned if variable names of TT1 and 
%   TT2 are the same; otherwise false is returned.
%
%   When called without an output, no error is thrown and the function silently
%   runs if variable names of TT1 and TT2 are consistent; otherwise an error 
%   is thrown to indicate inconsistency of variable names of TT1 and TT2.
%
% Inputs:
%   TT1 - First table or timetable tested for consistency. 
%
%   TT2 - Second table or timetable tested for consistency. 
%
%   Names - String vector or cell vector of character vectors ("cellstr")
%     with two elements indicating the names of TT1 and TT2. These names are
%     used to format more user-friendly error messages.
%
% Optional Output:
%   I - Logical indicator: true indicates that the variables of TT1 and TT2 
%     are consistent; false otherwise.
%
% Example:
%   Suppose TT1 and TT2 are timetables documented to users as "Y0" and "Y". 
%   This method is called as follows:
%
%   areTableVariableNamesConsistent(TT1,TT2,{'Y0' 'Y'})
%
% Note:
%   If both TT1 and TT2 are doubles/numerics, then no comparison is possible
%   and the indicator is defined to be true.
%

if isnumeric(TT1) || isnumeric(TT2)
   if nargout == 1
      varargout = {true};
   end
   return
end

%
% Ensure that tables/timetables have the same variable names.
%

I = isequal(TT1.Properties.VariableNames, TT2.Properties.VariableNames);
if nargout == 0
   if ~I
      error(message('econ:internal:econ:TableAndTimeTableUtilities:InconsistentTableVariableNames',names{1},names{2}))
   end
else
   varargout = {I};
end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = areTimeTableDateTimesConsistent(TT1,TT2,ConsistencyType,names)
% True if row dates/times of two timetables are consistent
%
% Syntax:
%   I = areTimeTableDateTimesConsistent(TT1,TT2,ConsistencyType,Names)
%   areTimeTableDateTimesConsistent(TT1,TT2,ConsistencyType,Names)
%
% Description:
%   Given regular timetables TT1 and TT2, determine if the row dates/times 
%   are "consistent" in one of two ways:
%
%    o Row times of TT1 projected into the future equal the row times of TT2
%      such that the concatenated row times of TT1 and TT2 represent a regular
%      sequence of date/time stamps. In other words, TT1 and TT2 are sequential
%      and occur in tandem, as if TT1 is the "pre-sample" data and TT2 the 
%      "in-sample" data of the same underlying time series Y(t) and such that
%      TT2 represents the future evolution or continuation of TT1. 
%
%    o The most recent row times of the longer of TT1 and TT2 are equal to
%      the row times of the shorter of TT1 and TT2 such that the most recent 
%      data in the longer series are concurrent with data in the shorter. 
%      For example, TT1 might represent a predictor series X(t) with 110 
%      observations and TT2 a response series Y(t) with 100 observations.
%      This function ensures that the last 100 row times of TT1 (X) coincide 
%      with the row times of TT2 (Y). 
%
%   When called with an output, true is returned if TT1 and TT2 are date/time
%   consistent; otherwise false is returned.
%
%   When called without an output, no error is thrown and the function silently
%   runs if TT1 and TT2 are consistent; otherwise an error is thrown to indicate
%   inconsistency of TT1 and TT2.
%
% Inputs:
%   TT1 - First timetable tested for consistency. 
%
%   TT2 - Second timetable tested for consistency. 
%
%   ConsistencyType - String or character vector indicating the method by 
%     which timetable row time consistency is determined:
%
%     o 'Sequential'  Row times of TT1 projected into the future equal the 
%           row times of TT2. When this option is selected TT1 must precede
%           TT2 in time.
%
%     o 'Concurrent'  The most recent row times of the longer of TT1 and TT2
%           are equal to the row times of the shorter of TT1 and TT2. This
%           option ensures row time consistency of the most recent observations
%           of TT1 and TT2.
%
%   Names - String vector or cell vector of character vectors ("cellstr")
%     with two elements indicating the names of TT1 and TT2. These names are
%     used to format more user-friendly error messages.
%
% Optional Output:
%   I - Logical indicator: true indicates that the row times of TT1 and TT2 
%     are consistent; false otherwise.
%
% Examples:
%   o The following regular timetables are sequentially row-time-consistent
%     because the 4 quarterly dates following the last date of TT1 (01-Apr-1960) 
%     appear as the first 4 quarterly dates in TT2 (01-Jul-1960, 01-Oct-1960, 
%     01-Jan-1961, and 01-Apr-1961):
%
%        Time           Y1           Y2           Y3   
%     ___________    _________    ________    ________
%     01-Jan-1960    -0.00557     0.030570     0.014354
%     01-Apr-1960     0.03297     0.042111     0.030412
%
%        Time           Y1           Y2           Y3    
%     ___________    ________     ________    _________
%     01-Jul-1960     0.03714     0.016360     0.031749
%     01-Oct-1960     0.09436     0.031930     0.024257
%     01-Jan-1961    -0.04359     0.021381    -0.002181
%     01-Apr-1961     0.02445     0.001921     0.044831
%
%   o The following regular timetables are concurrently row-time-consistent
%     because the last 4 quarterly dates of TT1 and are the same as the dates
%     of TT2:
%
%        Time           X1          X2          X3    
%     ___________    ________    ________    _________
%     01-Oct-1961     465.4       16.104       258.7
%     01-Jan-1962     468.9       16.276       259.8
%     01-Apr-1962     470.6       16.485       260.6 
%     01-Jul-1962     472.8       16.601       262.5 
%     01-Oct-1962     480.3       16.701       265.1 
%     01-Jan-1963     475.7       16.711       263.7 
%
%        Time           Y1          Y2          Y3    
%     ___________    ________    ________    _________
%     01-Apr-1962    0.021599    0.028270     0.011696
%     01-Jul-1962    0.012739    0.015558     0.017291
%     01-Oct-1962   -0.140180    0.013629     0.007590
%     01-Jan-1963    0.193580    0.013446     0.016870
%
% Notes:
%   o If neither TT1 nor TT2 is a timetable, then no comparison of row times 
%     is possible and the indicator is defined to be true.
%
%   o When TT1 is an empty timetable, sequential consistency of TT1 and TT2 
%     is ensured by assigning TT2's RowTimes to those of TT1. This is designed
%     to handle edge cases in which no pre-sample responses are needed, such 
%     as when simulating a VAR(0) model.
%

%
% Check dates/times if TT1 and TT2 are both timetables.
%

if ~istimetable(TT1) || ~istimetable(TT2)
   if nargout == 1
      varargout = {true};
   end
   return
end

if strcmpi(ConsistencyType, 'Sequential')

   if isempty(TT1.Properties.RowTimes)
%
%     Handle the edge case in which no RowTimes exist (e.g., simulating a
%     VAR(0) model in which no pre-sample responses are needed) separately
%     by simply assigning TT2's RowTimes to TT1's RowTimes. Since no pre-sample
%     data is needed, this forces sequential consistency and avoids an indexing 
%     error that would result if the assignment in the ELSE segment is executed.
%
      dates2 = TT2.Properties.RowTimes;
      dates1 = dates2;
   else
      if isnan(TT1.Properties.TimeStep)
%
%        Compute dates from TT2 into the past.
%
         dates1 = TT1.Properties.RowTimes;
         dates2 = TT2.Properties.RowTimes(1) - TT2.Properties.TimeStep * (size(TT1,1):-1:1)';
      else
%
%        Compute dates from TT1 into the future.
%
         dates1 = TT1.Properties.RowTimes(end) + TT1.Properties.TimeStep * (1:size(TT2,1))';
         dates2 = TT2.Properties.RowTimes;
      end
   end
elseif strcmpi(ConsistencyType, 'Concurrent')
%
%  Assign dates of TT1 & TT2 presumed to be concurrent.
%
   dates1 = TT1.Properties.RowTimes;
   dates2 = TT2.Properties.RowTimes;

   if size(TT1,1) >= size(TT2,1)
      dates1 = dates1(max(numel(dates1) - numel(dates2) + 1,1):end);
   else
      dates2 = dates2(max(numel(dates2) - numel(dates1) + 1,1):end);
   end
end

I = isequal(dates1, dates2);

if nargout == 0
   if ~I
      error(message('econ:internal:econ:TableAndTimeTableUtilities:InconsistentTimeTableDateTimes',names{1},names{2},lower(ConsistencyType)))
   end
else
   varargout = {I};
end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = isTimeTableRegular(TT,name)
% True if timetable is determined to be regular
%
% Syntax:
%   I = isTimeTableRegular(TT,Name)
%   isTimeTableRegular(TT,Name)
%
% Description:
%   Given a timetable TT, determine if TT satisfies all of the following:
%
%   (1) When TT has more than one observation (row), TT.Properties.TimeStep 
%       must be non-NaN (generally a calendar duration) so that TT is regular 
%       and past/future dates can be computed. Exceptions are made when TT 
%       is empty or has only a single observation, in which case TT is defined 
%       to be regular.
%   (2) TT data must be sorted in ascending order
%   (3) TT must have no duplicate dates/times
%   (4) TT.Properties.RowTimes and TT.Properties.StartTime must be datetimes
%
%   When called with an output, true is returned if all of the conditions are
%   satisfied; otherwise false is returned.
%
%   When called without an output, no error is thrown and the function silently
%   runs if all of the conditions are satisfied; otherwise an error is thrown.
%
% Inputs:
%   TT - Timetable tested for the conditions above.
%
%   Name - String scalar, character vector, or a scalar cell of character 
%     vector (1-by-1 "cellstr") indicating the name of TT in the caller's 
%     workspace. This name is used to a format more user-friendly error message.
%
% Optional Output:
%   I - Logical indicator: true indicates that TT satisfies the conditions
%     above; false otherwise.
%
% Notes:
%   o This utility may not be fully necessary after business calendar
%     awareness is added to datetime. Moreover, as of R2019a it is my 
%     understanding that g1846300 enhancement to timetable/isregular method 
%     would eliminate steps (1) and (3) in the list above; steps (2) and (4) 
%     would still need to be enforced.
%
%   o If TT is a conventional numeric array or a table (i.e., not a timetable),
%     then the indicator I is defined to be true.
%

if ~istimetable(TT)
   if nargout > 0
      varargout = {true};
   end
   return
end

%
% Get time dimension name so TT referencing is agnostic to user preferences.
%

Time = TT.Properties.DimensionNames{1};                  % TT default is 'Time'

%
% Test for the conditions above.
%

I = (size(TT,1) <= 1) || ~isnan(TT.Properties.TimeStep); % Is TimeStep non-NaN?
I = I && issorted(TT);                                   % and is TT sorted?
I = I && isequal(TT.(Time),unique(TT.(Time)));           % and are there no duplicate dates/times?
I = I && isdatetime(TT.Properties.StartTime);            % and is start time (and row times) a datetime?

if nargout == 0
   if ~I
      error(message('econ:internal:econ:TableAndTimeTableUtilities:IrregularTimeTable',char(name)))
   end
else
   varargout = {I};
end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = isTableDataNumericWithOutMissingData(TT,name)
% True if table/timetable data is numeric and has no missing data (NaNs)
%
% Syntax:
%   I = isTableDataNumericWithOutMissingData(TT,Name)
%   isTableDataNumericWithOutMissingData(TT,Name)
%
% Description:
%   Given table/timetable TT, determine if the underlying time series data is 
%   numeric (in most cases double) and contains no missing observations (NaNs).
%
%   When called with an output, true is returned if TT is numeric and without
%   missing data; otherwise false is returned.
%
%   When called without an output, no error is thrown and the function silently
%   runs if TT is numeric without missing data; otherwise an error is thrown.
%
% Inputs:
%   TT - Table or timetable tested for numeric data with no NaNs.
%
%   Name - String scalar, character vector, or a scalar cell of character 
%     vector (1-by-1 "cellstr") indicating the name of TT in the caller's 
%     workspace. This name is used to a format more user-friendly error message.
%
% Optional Output:
%   I - Logical indicator: true indicates that the data in TT is numeric and 
%   without missing data; false otherwise.
%
% Note:
%   If TT is a conventional numeric array (i.e., not a table or timetable),
%   then the indicator I is defined to be true and TT is not tested.
%

if isnumeric(TT)
   if nargout > 0
      varargout = {true};
   end
   return
end

I = isnumeric(TT.(TT.Properties.DimensionNames{2}));                % Is it all numeric?
if I
   I = I && all(all(~isnan(TT.(TT.Properties.DimensionNames{2})))); % And with no NaNs?
end

if nargout == 0
   if ~I
      error(message('econ:internal:econ:TableAndTimeTableUtilities:NonNumericOrMissingData',char(name)))
   end
else
   varargout = {I};
end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = isTableFormatValid(TT,name)
% True if table/timetable data is not nested inside another table/timetable
% and that each table/timetable variable has the same number of paths  
%
% Syntax:
%   I = isTableFormatValid(TT,Name)
%   isTableFormatValid(TT,Name)
%
% Description:
%   Given a table/timetable TT, determine if the underlying time series data 
%   is stored in a "flat" or "non-nested" format such that the data is not 
%   stored inside an inner table stored inside an outer table/timetable. Also, 
%   determine if each variable of TT has the same number of paths (columns).
%
%   Additional format validation checks may be added in the future.
%
%   When called with an output, true is returned if the format of TT is valid;
%   otherwise false is returned.
%
%   When called without an output, no error is thrown and the function silently
%   runs if the format of TT is valid; otherwise an error is thrown to indicate
%   that the format of TT is invalid.
%
% Inputs:
%   TT - Table or timetable tested for "flatness" or "non-nestedness" and 
%     variable path consistency.
%
%   Name - String scalar, character vector, or a scalar cell of character 
%     vector (1-by-1 "cellstr") indicating the name of TT in the caller's 
%     workspace. This name is used to a format more user-friendly error message.
%
% Optional Output:
%   I - Logical indicator: true indicates that the format of TT is valid; 
%     false otherwise.
%
% Examples:
%   Consider the following timetables:
%
%    (1) The data in the following TT is not nested and each variable has 3 
%        paths (columns), and so true is returned:
%
%          Time                  GDP                        COE          
%       ___________    _______________________    _______________________
%       31-Mar-2017    19065    19074    19062    10360    10348    10359
%       30-Jun-2017    19368    19248    19266    10447    10482    10409
%       30-Sep-2017    19488    19442    19452    10560    10569    10532
%       31-Dec-2017    19653    19615    19582    10587    10694    10545
%
%    (2) The data in the following TT is nested and false is returned:
%
%          Time                 Path1                      Path2         
%                       GDP      COE     PCEC      GDP      COE     PCEC 
%       ___________    _______________________    _______________________
%       31-Mar-2017    19065    10360    13190    19074    10348    13129
%       30-Jun-2017    19368    10447    13412    19248    10482    13232
%       30-Sep-2017    19488    10560    13473    19442    10569    13387
%       31-Dec-2017    19653    10587    13652    19615    10694    13549
%
%    (3) The data in the following TT is not nested but GDP has 2 paths and
%        COE has 3 paths, and so false is returned:
%
%          Time              GDP                   COE          
%       ___________    ______________    _______________________
%       31-Mar-2017    19065    19074    10360    10348    10359
%       30-Jun-2017    19368    19248    10447    10482    10409
%       30-Sep-2017    19488    19442    10560    10569    10532
%       31-Dec-2017    19653    19615    10587    10694    10545
%
% Note:
%   If TT is a conventional numeric array (i.e., not a table or timetable),
%   then the indicator I is defined to be true and TT non-nested.
%

if isnumeric(TT)
   if nargout > 0
      varargout = {true};
   end
   return
end

I = true;

%
% Ensure that table/timetable data is not nested inside another table/timetable.
%

for i = 1:size(TT,2)
    I = I && isnumeric(TT.(TT.Properties.VariableNames{i}));
    if ~I 
       break   % as soon as it becomes false
    end
end

%
% The variables are numeric, so now ensure that they have the same number of paths.
%

if I
   I = I && isnumeric(TT.(TT.Properties.DimensionNames{2}));
   for i = 2:size(TT,2)
       I = I && ( size(TT(:,1).(TT.Properties.DimensionNames{2}),2) == size(TT(:,i).(TT.Properties.DimensionNames{2}),2) );
       if ~I 
          break   % as soon as it becomes false
       end
   end
end

if nargout == 0
   if ~I
      error(message('econ:internal:econ:TableAndTimeTableUtilities:InvalidTableFormat',char(name)))
   end
else
   varargout = {I};
end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = isTableDataSinglePath(TT,name)
% True if table/timetable stores a single path of time series data
%
% Syntax:
%   I = isTableDataSinglePath(TT,Name)
%   isTableDataSinglePath(TT,Name)
%
% Description:
%   For a table or timetable TT, determine if all variables in TT store a 
%   single path of time series data.
%
%   When called with an output, true is returned if all variable in TT store
%   a single column vector of numeric data; otherwise false is returned.
%
%   When called without an output, no error is thrown and the function silently
%   runs if all variables in TT store a single path; otherwise an error is 
%   thrown.
%
% Inputs:
%   TT - Table or timetable tested to determine if all variables store a single
%     path of time series data.
%
%   Name - String scalar, character vector, or a scalar cell of character 
%     vector (1-by-1 "cellstr") indicating the name of TT in the caller's 
%     workspace. This name is used to a format more user-friendly error message.
%
% Optional Output:
%   I - Logical indicator: true indicates that all variables in TT store a 
%     single path of time series data; false otherwise.
%
% Examples:
%   Consider the following timetables:
%
%    (1) The following TT has 2 variables and each variable has a single 
%        path (column) of numeric data, and so true is returned:
%
%           Time         GDP      COE          
%        ___________    _____    _____
%        31-Mar-2017    19065    10360
%        30-Jun-2017    19368    10447
%        30-Sep-2017    19488    10560
%        31-Dec-2017    19653    10587
%
%    (2) The following TT has 2 variables and each variable has 3 paths of 
%        numeric data, and so false is returned:
%
%           Time                  GDP                        COE          
%        ___________    _______________________    _______________________
%        31-Mar-2017    19065    19074    19062    10360    10348    10359
%        30-Jun-2017    19368    19248    19266    10447    10482    10409
%        30-Sep-2017    19488    19442    19452    10560    10569    10532
%        31-Dec-2017    19653    19615    19582    10587    10694    10545
%
%    (2) The following TT has 2 variables and although the first variable has
%        a single path of numeric data, the second variable (COE) has 2 paths,
%        and so false is returned:
%
%           Time         GDP           COE          
%        ___________    _____    ______________
%        31-Mar-2017    19065    10360    10348
%        30-Jun-2017    19368    10447    10482
%        30-Sep-2017    19488    10560    10569
%        31-Dec-2017    19653    10587    10694
%
% Note:
%   If TT is a conventional numeric array (i.e., not a table or timetable),
%   then the indicator I is defined to be true and TT is not tested to ensure
%   that all variables store a single path. 
%
%   To test whether numeric time series TT represents a single path requires 
%   an input flag to indicate whether TT represents a univariate or multivariate
%   time series. For numeric, univariate TT, a single path requires that 
%   size(TT,2) = 1, while multivariate TT requires size(TT,3) = 1.
%

if isnumeric(TT)
   if nargout > 0
      varargout = {true};
   end
   return
end

for i = 1:size(TT,2)
    I = (size(TT(:,i).(TT.Properties.DimensionNames{2}),2) == 1);
    if ~I 
       break   % as soon as it becomes false
    end
end

if nargout == 0
   if ~I
      error(message('econ:internal:econ:TableAndTimeTableUtilities:NonSinglePathData',char(name)))
   end
else
   varargout = {I};
end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function validateTablesAndTimeTables(varargin)
% Ensure tables/timetables are not nested, contain numeric data with no missing
% observations, and future dates can be computed
%
% Syntax:
%   validateTablesAndTimeTables(TT1,TT2,...,TTN,Names)
%
% Description:
%   Validate an arbitrary list of tables/timetables to ensure that all 
%   tables/timetables are not nested inside other tables/timetables, contain 
%   numeric data with no missing observations (NaNs), and that all timetables
%   are determined to be regular.
%
%   This is a convenience wrapper designed to minimize code changes to toolbox 
%   functions around the following fundamental utilities:
%
%     o isTimeTableRegular
%     o isTableFormatValid
%     o isTableDataNumericWithOutMissingData
%
%   This utility has no outputs, but errors are thrown if any of the inputs
%   fail validation.
%
% Inputs:
%   TT1, TT2, ... TTN - List of tables and timetables to validate.
%
%   Names - String vector or cell vector of character vectors ("cellstr")
%     of the published names of TT1, TT2, ..., TTN in the caller's workspace 
%     (i.e., the published names of the variables in the calling function 
%     associated with TT1, TT2, ..., TTN). 
%
% Note:
%   Any of TT1, TT2, ..., TTN may be a conventional numeric array (i.e., not 
%   a table or timetable), in which case no validation is performed.
%

names    = string(varargin{end});          % The last element stores variable names
varargin = varargin(1:end-1);              % Retain only time series inputs TT1, TT2, ..., TTN

for i = 1:numel(varargin)
    try 
      internal.econ.TableAndTimeTableUtilities.isTableFormatValid(varargin{i},names(i))
    catch exception
      exception.throwAsCaller()
    end
    try
      internal.econ.TableAndTimeTableUtilities.isTimeTableRegular(varargin{i},names(i))
    catch exception
      exception.throwAsCaller()
    end
    try
      internal.econ.TableAndTimeTableUtilities.isTableDataNumericWithOutMissingData(varargin{i},names(i))
    catch exception
      exception.throwAsCaller()
    end    
end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function ensureTimeSeriesTypeConsistency(varargin)
% Ensure time series inputs are the same data type
%
% Syntax:
%   ensureTimeSeriesTypeConsistency(TT1, TT2, ..., TTN, Names, UsingDefaults)
%
% Description:
%   Enforce type consistency between time series inputs, ensuring that all 
%   are either doubles/numerics, or all tables, or all timetables.
%
%   This utility has no outputs, but an error is thrown if the data type of 
%   time series TT1, TT2, ..., TTN differ.
%
% Inputs:
%   TT1, TT2, ..., TTN - List of time series tested for type consistency.
%
%   Names - String vector or cell vector of character vectors ("cellstr")
%     of the published names of TT1, TT2, ..., TTN in the caller's workspace 
%     (i.e., the published names of the variables in the calling function 
%     associated with TT1, TT2, ..., TTN). 
%
%   UsingDefaults - String vector or cell vector of character vectors ("cellstr")
%     of names of optional input N-V pairs published in the caller's workspace 
%     using their default values (i.e., published N-V pairs not explicitly 
%     specified by the user). Variable names in UsingDefaults are not tested 
%     for type consistency because they are unspecified by the user.
%
% Note:
%   o It is typical that the calling function/method uses "inputParser" to
%     validate inputs, in which case UsingDefaults = parser.UsingDefaults.
%
% Example:
%   Suppose function FOO accepts optional time series inputs 'Y0', 'Y', and
%   'X', and that the user calls it with only 'Y' and 'X':
%
%   [...] = FOO(..., 'Y', Y, 'X', X, ...);
%
%   For input parser P, the UsingDefaults property will contain 'Y0', and so
%   the following tests for data type consistency between 'Y' and 'X' and
%   excludes 'Y0':
%
%   ensureTimeSeriesTypeConsistency(Y0, Y, X, {'Y0' 'Y' 'X'}, P.UsingDefaults)
%

usingDefaults = varargin{end};                  % The last element stores default information
names         = varargin{end-1};                % The 2nd-to-last element stores variable names
varargin      = varargin(1:end-2);              % Retain only time series inputs TT1, TT2, ..., TTN
type          = repmat("", numel(varargin), 1); % Pre-allocate data type vector

if numel(type) <= 1
   return      % There must be at least 2 time series to compare
end

iCount = 0;    % # of time series to compare

for i = 1:numel(varargin)
    iMatch  = strcmp(names(i), usingDefaults);
    if ~any(iMatch)
       iCount       = iCount + 1;
       type(iCount) = string(class(varargin{i}));
    end
end

if numel(unique(type(1:iCount))) > 1
   error(message('econ:internal:econ:TableAndTimeTableUtilities:InconsistentTimeSeriesDataTypes'))
end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function [Y,MetaData] = perVariableTable2PerPathMatrix(TT)
% Convert per variable table/timetable to per path 2-D univariate matrix
%
% Syntax:
%
%	[Y,MetaData] = perVariableTable2PerPathMatrix(TT)
%
% Description:
%
%   Given a table or timetable TT in "per variable" format, extract the data
%   and return a 2-D matrix Y in "per path" format in which each column 
%   represents a distinct path of the underlying univariate time series.
%
% Input:
%   TT - Table or timetable whose data is stored in per variable format. 
% 
% Outputs:
%   Y - 2-D matrix whose data is stored in per path format.
%
%   MetaData - Metadata property information of TT. MetaData is a 
%     'matlab.tabular.TableProperties' or 'matlab.tabular.TimetableProperties' 
%     object when TT is a table or timetable, respectively. 
%
% Example:
%   Suppose a timetable TT stores 3 paths of GDP over 4 quarters in per 
%   variable format:
%
%        Time                  GDP          
%     ___________    _______________________
%     31-Mar-2017    19067    19218    19034
%     30-Jun-2017    19353    19619    19294
%     30-Sep-2017    19452    19846    19560
%     31-Dec-2017    19654    20233    19871
%
%   Then the corresponding matrix Y in per path format is:
%
%     Y =
%         19067  19218  19034
%         19353  19619  19294
%         19452  19846  19560
%         19654  20233  19871
%
% Note:
%   If TT is a conventional numeric array (i.e., not a table or timetable),
%   then the input is simply passed through such that Y = TT.
%

if isnumeric(TT)
   Y        = TT;
   MetaData = [];
   return
end

Y = TT{:,:};

if nargout > 1
   MetaData = TT.Properties;
end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function TT = perPathMatrix2PerVariableTable(Y, MetaData, MetaDataPropertiesToExclude, DateTimeType)    
% Convert per path 2-D univariate matrix to per variable table/timetable
%
% Syntax:
%
%	TT = perPathMatrix2PerVariableTable(Y,MetaData,PropertiesToExclude,DateTimeType)
%
% Description:
%
%   Given a 2-D univariate time series matrix Y in "per path" format, extract
%   the data and return a table or timetable TT in "per variable" format. As
%   part of the conversion, certain metadata property information is also 
%   copied to TT (i.e., this utility also calls "copyMetaData2Table"; see 
%   notes below for details).
%
% Inputs:
%   Y - 2-D matrix whose data is stored in per path format.
%
%   MetaData - Metadata property information used to determine whether TT is
%     returned as a table or timetable. MetaData is a 'matlab.tabular.TableProperties'
%     or 'matlab.tabular.TimetableProperties' object and TT is then a table 
%     or timetable, respectively. MetaData also stores variables names used
%     to create both tables and timetables, as well as time/date stamp
%     information used to populate row times of timetables. When Y is a table 
%     or timetable, MetaData is ignored.
%
%   PropertiesToExclude - List of metadata properties excluded from the copy 
%     operation, specified as a string vector or a cell vector of character 
%     vectors ("cellstr"). TT properties in this list are not over-written 
%     by properties in MetaData. See notes below for additional details.
%
%   DateTimeType - String or character vector indicating the method by which 
%     timetable time/date stamps are determined. For 2-D array Y with N
%     observations:
%
%     o 'First'  Row times are taken directly from the first N row times of
%                MetaData:
%
%                MetaData.RowTimes(1:N)
%
%     o 'Last'   Row times are taken directly from the last N row times of 
%                MetaData:
%
%                MetaData.RowTimes(numel(MetaData.RowTimes) - N + 1:end)
%
%     o 'Next'   Row times begin at the next datetime immediately following 
%                last datetime of MetaData (i.e., datetimes in the future):
%
%                MetaData.RowTimes(end) + MetaData.TimeStep * (1:N)
%
% Outputs:
%   TT - Table or timetable whose data is stored in per variable format. 
%
% Example:
%   Suppose a time series matrix Y stores 3 paths of GDP over 4 quarters in 
%   per path format:
%
%     Y =
%         19067  19218  19034
%         19353  19619  19294
%         19452  19846  19560
%         19654  20233  19871
%
%   The corresponding timetable TT in per variable format is, e.g.,:
%
%        Time                  GDP          
%     ___________    _______________________
%     31-Mar-2017    19067    19218    19034
%     30-Jun-2017    19353    19619    19294
%     30-Sep-2017    19452    19846    19560
%     31-Dec-2017    19654    20233    19871
%
% Notes:
%   o If Y is a table or timetable (i.e., not a conventional numeric array),
%     then the input is simply passed through such that TT = Y. Also, an empty
%     MetaData indicates that Y was originally specified as a numeric array 
%     and that no conversion should take place, and again the input is simply 
%     passed through such that TT = Y.
%
%   o Metadata Property Copy Details:
%     - For both tables and timetables, metadata property "VariableNames" is 
%       always copied to TT such that variable names are preserved.
%
%     - For timetables, date/time stamps of TT are derived from MetaData and 
%       the input choice of DateTimeType (see discussion above). Since TT
%       is a regular timetable, "TimeStep" (i.e., MetaData.TimeStep) is always
%       copied to TT, but "StartTime", "RowTimes", and "SampleRate" are derived.
%
%   o When MetaData.RowTimes is empty, TT.RowTimes begins at the start time 
%     of MetaData (i.e., MetaData.StartTime). This is designed to handle 
%     edge cases in which no pre-sample responses are needed, such as when 
%     simulating a VAR(0) model.

if isempty(MetaData) || ~isnumeric(Y)
   TT = Y;
   return
end

nDates = size(Y,1);

if isa(MetaData,'matlab.tabular.TableProperties')    % Is Y a Table?
   TT = table('Size', [nDates 1], 'VariableTypes', {class(Y)}, ...
              'VariableNames', MetaData.VariableNames);
else                                                 % Y is a Timetable
   if strcmpi(DateTimeType, 'Next')
%
%     Compute dates into the future (i.e., the next dates).
%
      if isempty(MetaData.RowTimes)
%
%        In the edge case in which no RowTimes exist (e.g., simulating a
%        VAR(0) model in which no pre-sample responses are needed), begin
%        TT dates at the StartTime found in MetaData. This is a convention
%        which, lacking any additional information, makes as much sense as
%        anything.
%
         startTime = MetaData.StartTime;
      else
         startTime = MetaData.RowTimes(end) + MetaData.TimeStep;
      end
   elseif strcmpi(DateTimeType, 'First')
%
%     Assign dates from beginning (i.e., the first dates).
%
      startTime = MetaData.StartTime;
   elseif strcmpi(DateTimeType, 'Last')
%
%     Assign dates from end (i.e., the last dates).
%
      if nDates <= 1
%
%        Account for the edge case in which we want an output TT with only 
%        one row but whose MetaData.TimeStep = NaN.
%
         startTime = MetaData.StartTime;
      else
         startTime = MetaData.StartTime + MetaData.TimeStep * (numel(MetaData.RowTimes) - nDates);
      end
   end

   TT = timetable('Size', [nDates 1], 'VariableTypes'       , {class(Y)} , ...
                  'VariableNames'   , MetaData.VariableNames, 'StartTime', startTime, ...
                  'TimeStep'        , MetaData.TimeStep);
end

TT.(MetaData.VariableNames{1}) = Y;

%
% Copy remaining metadata information. 
%
% Since "VariableNames", "RowTimes", "StartTime", "SampleRate", and "TimeStep" 
% have already been determined, they are excluded from the direct copy operation.
%

MetaDataPropertiesToExclude = unique(cat(1,MetaDataPropertiesToExclude(:), ...
                                    {'VariableNames' 'RowTimes' 'StartTime' 'SampleRate' 'TimeStep'}'));
TT = internal.econ.TableAndTimeTableUtilities.copyMetaData2Table(TT,MetaData,MetaDataPropertiesToExclude);

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function [Y,MetaData] = perVariableTable2PerPathArray(TT)
% Convert per variable table/timetable to per path 3-D multivariate array
%
% Syntax:
%
%	Y = perVariableTable2PerPathArray(TT)
%
% Description:
%
%   Given a table or timetable TT in "per variable" format, extract the data
%   and return a 3-D array Y in "per path" format.
%
% Inputs:
%   TT - Table or timetable whose data is stored in per variable format. 
% 
% Outputs:
%   Y - 3-D array whose data is stored in per path format.
%
%   MetaData - Metadata property information of TT. MetaData is a 
%     'matlab.tabular.TableProperties' or 'matlab.tabular.TimetableProperties' 
%     object when TT is a table or timetable, respectively. 
%
% Example:
%   Suppose a timetable TT stores 3 paths of 2 variables over 4 quarters in
%   per variable format:
%
%        Time                  GDP                      GDPDEF         
%     ___________    _______________________    _______________________
%     31-Mar-2017    19067    19218    19034    10384    10432    10388
%     30-Jun-2017    19353    19619    19294    10470    10638    10454
%     30-Sep-2017    19452    19846    19560    10584    10835    10637
%     31-Dec-2017    19654    20233    19871    10657    10965    10830
%
%   Then the corresponding 3-D array Y in per path format is:
%
%     Y(:,:,1) =
%         19067        10384
%         19353        10470
%         19452        10584
%         19654        10657
%     Y(:,:,2) =
%         19218        10432
%         19619        10638
%         19846        10835
%         20233        10965
%     Y(:,:,3) =
%         19034        10388
%         19294        10454
%         19560        10637
%         19871        10830
%
% Note:
%   If TT is a conventional numeric array (i.e., not a table or timetable),
%   then the input is simply passed through such that Y = TT and MetaData is 
%   empty.
%

if isnumeric(TT)
   Y        = TT;
   MetaData = [];
   return
end

[nDates,nVariables] = size(TT);
nPaths              = size(TT{:,1},2);
Y                   = zeros(nDates,nVariables,nPaths);

for i = 1:nPaths
    Y(:,:,i) = TT{:,:}(:,i:nPaths:end);
end

if nargout > 1
   MetaData = TT.Properties;
end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function TT = perPathArray2PerVariableTable(Y, MetaData, MetaDataPropertiesToExclude, DateTimeType)
% Convert per path 3-D multivariate array to per variable table/timetable
%
% Syntax:
%
%	TT = perPathArray2PerVariableTable(Y,MetaData,PropertiesToExclude,DateTimeType)
%
% Description:
%
%   Given a 3-D multivariate time series array Y in "per path" format, extract
%   the data and return a table or timetable TT in "per variable" format. As
%   part of the conversion, certain metadata property information is also 
%   copied to TT (i.e., this utility also calls "copyMetaData2Table"; see 
%   notes below for details).
%
% Inputs:
%   Y - 3-D multivariate array whose data is stored in per path format.
%
%   MetaData - Metadata property information used to determine whether TT is
%     returned as a table or timetable. MetaData is a 'matlab.tabular.TableProperties'
%     or 'matlab.tabular.TimetableProperties' object and TT is then a table 
%     or timetable, respectively. MetaData also stores variables names used
%     to create both tables and timetables, as well as time/date stamp
%     information used to populate row times of timetables. When Y is a table 
%     or timetable, MetaData is ignored.
%
%   PropertiesToExclude - List of metadata properties excluded from the copy 
%     operation, specified as a string vector or a cell vector of character 
%     vectors ("cellstr"). TT properties in this list are not over-written 
%     by properties in MetaData. See notes below for additional details.
%
%   DateTimeType - String or character vector indicating the method by which 
%     timetable TT row time stamps are determined. For 3-D array Y with N
%     observations:
%
%     o 'First'  Row times are taken directly from the first N row times of
%                MetaData:
%
%                MetaData.RowTimes(1:N)
%
%     o 'Last'   Row times are taken directly from the last N row times of 
%                MetaData:
%
%                MetaData.RowTimes(numel(MetaData.RowTimes) - N + 1:end)
%
%     o 'Next'   Row times begin at the next datetime immediately following 
%                last datetime of MetaData (i.e., datetimes in the future):
%
%                MetaData.RowTimes(end) + MetaData.TimeStep * (1:N)
%
% Outputs:
%   TT - Table or timetable whose data is stored in per variable format. 
%
% Example:
%   Suppose a 3-D array Y stores 3 paths of 2 variables over 4 quarters in 
%   per path format:
%
%     Y(:,:,1) =
%         19067        10384
%         19353        10470
%         19452        10584
%         19654        10657
%     Y(:,:,2) =
%         19218        10432
%         19619        10638
%         19846        10835
%         20233        10965
%     Y(:,:,3) =
%         19034        10388
%         19294        10454
%         19560        10637
%         19871        10830
%
%   Then the corresponding timetable TT in per variable format is, e.g.:
%
%        Time                  GDP                      GDPDEF         
%     ___________    _______________________    _______________________
%     31-Mar-2017    19067    19218    19034    10384    10432    10388
%     30-Jun-2017    19353    19619    19294    10470    10638    10454
%     30-Sep-2017    19452    19846    19560    10584    10835    10637
%     31-Dec-2017    19654    20233    19871    10657    10965    10830
%
% Notes:
%   o If Y is a table or timetable (i.e., not a conventional numeric array),
%     then the input is simply passed through such that TT = Y. Also, an empty
%     MetaData indicates that Y was originally specified as a numeric array 
%     and that no conversion should take place, and again the input is simply 
%     passed through such that TT = Y.
%
%   o Metadata Property Copy Details:
%     - For both tables and timetables, metadata property "VariableNames" is 
%       always copied to TT such that variable names are preserved.
%
%     - For timetables, date/time stamps of TT are derived from MetaData and 
%       the input choice of DateTimeType (see discussion above). Since TT
%       is a regular timetable, "TimeStep" (i.e., MetaData.TimeStep) is always
%       copied to TT, but "StartTime", "RowTimes", and "SampleRate" are derived.
%
%   o When MetaData.RowTimes is empty, the TT.RowTimes begins at the start
%     time of MetaData (i.e., MetaData.StartTime). This is designed to handle 
%     edge cases in which no pre-sample responses are needed, such as when 
%     simulating a VAR(0) model.

if isempty(MetaData) || ~isnumeric(Y)
   TT = Y;
   return
end

[nDates,nVariables,~] = size(Y);
types                 = repmat({class(Y)},1,nVariables);

if isa(MetaData,'matlab.tabular.TableProperties')    % Is Y a Table?
   TT = table('Size', [nDates nVariables], 'VariableTypes', types, ...
              'VariableNames', MetaData.VariableNames);
else                                                 % Y is a Timetable
   if strcmpi(DateTimeType, 'Next')
%
%     Compute dates into the future (i.e., the next dates).
%
      if isempty(MetaData.RowTimes)
%
%        In the edge case in which no RowTimes exist (e.g., simulating a
%        VAR(0) model in which no pre-sample responses are needed), begin
%        TT dates at the StartTime found in MetaData. This is a convention
%        which, lacking any additional information, makes as much sense as
%        anything.
%
         startTime = MetaData.StartTime;
      else
         startTime = MetaData.RowTimes(end) + MetaData.TimeStep;
      end
   elseif strcmpi(DateTimeType, 'First')
%
%     Assign dates from beginning (i.e., the first dates).
%
      startTime = MetaData.StartTime;
   elseif strcmpi(DateTimeType, 'Last')
%
%     Assign dates from end (i.e., the last dates).
%
      if nDates <= 1
%
%        Account for the edge case in which we want an output TT with only 
%        one row but whose MetaData.TimeStep = NaN.
%
         startTime = MetaData.StartTime;
      else
         startTime = MetaData.StartTime + MetaData.TimeStep * (numel(MetaData.RowTimes) - nDates);
      end
   end

   TT = timetable('Size'         , [nDates nVariables]   , 'VariableTypes', types    , ...
                  'VariableNames', MetaData.VariableNames, 'StartTime'    , startTime, ...
                  'TimeStep'     , MetaData.TimeStep);
end

Y = permute(Y,[1 3 2]);   % Convert from per path to per variable format

for i = 1:nVariables
    TT.(MetaData.VariableNames{i}) = Y(:,:,i);
end

%
% Copy remaining metadata information.
%
% Since "VariableNames", "RowTimes", "StartTime", "SampleRate", and "TimeStep" 
% have already been determined, they are excluded from the direct copy operation.
%

MetaDataPropertiesToExclude = unique(cat(1,MetaDataPropertiesToExclude(:), ...
                                    {'VariableNames' 'RowTimes' 'StartTime' 'SampleRate' 'TimeStep'}'));
TT = internal.econ.TableAndTimeTableUtilities.copyMetaData2Table(TT,MetaData,MetaDataPropertiesToExclude);

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function TT = copyMetaData2Table(TT,MetaData,MetaDataPropertiesToExclude)
% Update table/timetable metadata
%
% Syntax:
%   TT = copyMetaData2Table(TT,MetaData,PropertiesToExclude)
%
% Description:
%   Given an input table/timetable TT, update select metadata properties of 
%   TT with those of MetaData (i.e., copy select properties from MetaData to 
%   TT, over-writing the metadata in TT).
%
% Input:
%   TT - Table or timetable whose metadata is updated.
%
%   MetaData - Metadata property information used to update TT. MetaData is 
%     a 'matlab.tabular.TableProperties' or 'matlab.tabular.TimetableProperties' 
%     object. When TT is a not table or timetable, MetaData is ignored.
%
%   PropertiesToExclude - List of metadata properties excluded from the copy 
%     operation, specified as a string vector or a cell vector of character 
%     vectors ("cellstr"). TT properties in this list are not over-written.
%
% Output:
%   TT - Updated table or timetable whose metadata has been copied from 
%        MetaData (i.e., over-written by MetaData).
%
% Example:
%   Exclude 'VariableUnits' metadata property as follows:
%
%   TT = copyMetaData2Table(TT, MetaData, {'VariableUnits'})
%
% Note:
%   If Y is a conventional numeric array (i.e., not a table or timetable),
%   then the input TT is unchanged.
%

if isempty(MetaData) || isnumeric(TT)
   return
end

%
% Extract all metadata properties of the input table/timetable EXCEPT those
% specifically excluded.
%

fields = setdiff(fieldnames(MetaData), MetaDataPropertiesToExclude);

%
% Assign all remaining metadata properties of the input table/timetable to 
% the output table/timetable.
%

for i = 1:numel(fields)
    TT.Properties.(fields{i}) = MetaData.(fields{i});
end

end

end % Methods (Static)

end % Class definition