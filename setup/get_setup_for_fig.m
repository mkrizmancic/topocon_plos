function [A, l2_refs, signal_refs, params, doPlot] = get_setup_for_fig(ID)
%GET_SETUP_FOR_FIG Get the initial topology and experimemnt setup to recreate figure in the paper.

l2_refs = [];
signal_refs = {};
params = struct();
doPlot = struct();

switch ID
    case 3
        n = 6;
        A = diag(ones(1, n-1),1) + diag(ones(1, n-1),-1);
        l2_refs = [ 2, 150, 400,  900, 1500;
                   -1,   1, 2.3, 0.26, 0.26];

        params.steps = 1500;

        doPlot.Lambda = true;
        doPlot.subSigma = true;

    case 4
        A = [0 1 1 1 0 0;
             1 0 0 0 1 1;
             1 0 0 0 0 1;
             1 0 0 0 0 0;
             0 1 0 0 0 0;
             0 1 1 0 0 0];
        params.steps = 1;
        doPlot.Graph = true;

    case 5
        A = [0 1 1 1 0 0;
             1 0 0 0 1 1;
             1 0 0 0 0 1;
             1 0 0 0 0 0;
             0 1 0 0 0 0;
             0 1 1 0 0 0];
        l2_refs = [ 2,  100;       % step
                   -1, 0.65];      % value
        signal_refs = {150, {1, 4, 0.9;
                             4, 1, 0.9};
                       350, {1, 4, 0.7;
                             4, 1, 0.7}};

        params.steps = 700;

        doPlot.Lambda = true;
        doPlot.K_l2 = true;
        doPlot.subLink = true;

    case 6
        A = [0 1 1 1 0 0;
             1 0 0 0 1 1;
             1 0 0 0 0 1;
             1 0 0 0 0 0;
             0 1 0 0 0 0;
             0 1 1 0 0 0];
        l2_refs = [ 2, 100, 300, 500;       % step
                   -1,   2,   3, 1.6];      % value

        params.steps = 800;

        doPlot.Lambda = true;

    case 7
        A = [0 1 1 1 0 0;
             1 0 0 0 1 1;
             1 0 0 0 0 1;
             1 0 0 0 0 0;
             0 1 0 0 0 0;
             0 1 1 0 0 0];
        l2_refs = [ 2, 100, 300, 500;       % step
                   -1,   2,   3, 1.6];      % value

        params.steps = 800;

        doPlot.Snapshots = true;

    case 8
        n = 6;
        A = [0 1 1 1 0 0;
             1 0 0 0 1 1;
             1 0 0 0 0 1;
             1 0 0 0 0 0;
             0 1 0 0 0 0;
             0 1 1 0 0 0];
        l2_refs = [ 2, 100, 800;       % step
                   -1,   2,   3];      % value
        signal_refs = {350, {3, 1:n, 0};
                       600, {3, 1:n, 1}};

        params.steps = 1000;

        doPlot.Lambda = true;

    case 9
        A = [0 1 1 1 1 0;
             1 0 0 0 1 1;
             1 0 0 1 0 1;
             1 0 1 0 1 0;
             1 1 0 1 0 1;
             0 1 1 0 1 0];
        l2_refs = [ 2, 100, 300;
                   -1,   3,   4];

        params.steps = 800;
        params.remove_zero_links = true;

        doPlot.Lambda = true;
        doPlot.avgMatrix = true;

    case 10
        A = [0 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
             1 0 1 1 0 0 0 0 0 0 1 0 0 0 1;
             1 1 0 1 1 0 1 0 0 1 0 0 0 0 0;
             0 1 1 0 0 1 0 0 1 0 0 0 0 0 0;
             0 0 1 0 0 1 1 0 0 0 0 0 0 1 0;
             0 0 0 1 1 0 1 1 1 0 0 1 0 0 0;
             0 0 1 0 1 1 0 1 1 0 0 0 1 1 0;
             0 0 0 0 0 1 1 0 0 0 0 0 0 0 0;
             0 0 0 1 0 1 1 0 0 0 0 0 0 0 0;
             0 0 1 0 0 0 0 0 0 0 0 0 1 0 0;
             0 1 0 0 0 0 0 0 0 0 0 0 1 1 0;
             0 0 0 0 0 1 0 0 0 0 0 0 1 1 0;
             0 0 0 0 0 0 1 0 0 1 1 1 0 0 1;
             0 0 0 0 1 0 1 0 0 0 1 1 0 0 1;
             0 1 0 0 0 0 0 0 0 0 0 0 1 1 0];
        params.steps = 1;
        doPlot.Graph = true;

    case 11
        n = 15;
        A = [0 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
             1 0 1 1 0 0 0 0 0 0 1 0 0 0 1;
             1 1 0 1 1 0 1 0 0 1 0 0 0 0 0;
             0 1 1 0 0 1 0 0 1 0 0 0 0 0 0;
             0 0 1 0 0 1 1 0 0 0 0 0 0 1 0;
             0 0 0 1 1 0 1 1 1 0 0 1 0 0 0;
             0 0 1 0 1 1 0 1 1 0 0 0 1 1 0;
             0 0 0 0 0 1 1 0 0 0 0 0 0 0 0;
             0 0 0 1 0 1 1 0 0 0 0 0 0 0 0;
             0 0 1 0 0 0 0 0 0 0 0 0 1 0 0;
             0 1 0 0 0 0 0 0 0 0 0 0 1 1 0;
             0 0 0 0 0 1 0 0 0 0 0 0 1 1 0;
             0 0 0 0 0 0 1 0 0 1 1 1 0 0 1;
             0 0 0 0 1 0 1 0 0 0 1 1 0 0 1;
             0 1 0 0 0 0 0 0 0 0 0 0 1 1 0];
        l2_refs = [ 2, 200, 1000;
                   -1, 1.4, 1.2];
        signal_refs = {350, {3, 1:n, 0;
                             9, 1:n, 0;
                             12, 1:n, 0};
                       1400, {12, 1:n, 1}};

        params.steps = 2000;
        params.norms_limit = 0.05;

        doPlot.Lambda = true;

end
end