%% Allows saving in parallel
% This helper function allows saving within a `parfor` loop.
% The name of the saved structure is forced to be 'LLE_save_obj' to ensure consistency.
function parsave(fname, save_obj)
  save(fname, 'save_obj')
end