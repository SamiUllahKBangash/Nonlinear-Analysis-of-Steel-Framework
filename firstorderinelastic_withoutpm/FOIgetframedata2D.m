function [nodes, elems, bcs, loads, E, I, A, Mp] = FOIgetframedata2D(filename)
S=fileread(filename);
np = regexp(S, '^Nodes:', 'once', 'lineanchors');
ep = regexp(S, '^Elements:', 'once', 'lineanchors');
bcp = regexp(S, '^bcs:', 'once', 'lineanchors');
nlp = regexp(S, '^Nodal loads:', 'once', 'lineanchors');
nodes = cell2mat( textscan(S(np : ep-1), '%f %f', 'HeaderLines', 1) );
elements = cell2mat( textscan(S(ep:bcp-1), '%f %f %f %f %f %f', 'HeaderLines', 1));
bcs = cell2mat( textscan(S(bcp:nlp-1), '%f %f %f', 'HeaderLines', 1));
loads = cell2mat( textscan(S(nlp:end), '%f %f %f', 'HeaderLines', 1));
elems = elements(:,1:2);
E = elements(:,3);
A = elements(:,4);
I = elements(:,5);
Mp = elements(:,6);
end




