function files = remove_first2rows(files)

rowToDelete = [1 2]; % or whatever....
files(rowToDelete, :) = [];
