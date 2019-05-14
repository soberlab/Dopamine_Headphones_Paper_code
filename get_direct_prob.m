function [ p_test, p_joint_matrix ] = get_direct_prob( sample1, sample2 )
%get_direct_prob Returns the direct probability of items from sample2 being
%greater than or equal to those from sample1. 
%   Sample1 and Sample2 are two bootstrapped samples and this function
%   directly computes the probability of items from sample 2 being greater
%   than or equal to those from sample1. Since the bootstrapped samples are
%   themselves posterior distributions, this is a way of computing a
%   Bayesian probability. The joint matrix can also be returned to compute
%   directly upon.

%First, we want to find the ranges over which we need to build the
%probability distributions.
joint_low_val = min([min(sample1) min(sample2)]);
joint_high_val = max([max(sample1) max(sample2)]);

%Now set up a matrix with axis between these extreme values for which we
%will calculate the probabilities.

p_joint_matrix = zeros(100,100);
p_axis = linspace(joint_low_val,joint_high_val,100);
edge_shift = (p_axis(2) - p_axis(1))/2;
p_axis_edges = p_axis - edge_shift;
p_axis_edges = [p_axis_edges (joint_high_val + edge_shift)];

%Calculate probabilities using histcounts for edges.

p_sample1 = histcounts(sample1,p_axis_edges)/length(sample1);
p_sample2 = histcounts(sample2,p_axis_edges)/length(sample2);

%Now, calculate the joint probability matrix:

for i = 1:size(p_joint_matrix,1)
    for j = 1:size(p_joint_matrix,2)
        p_joint_matrix(i,j) = p_sample1(i)*p_sample2(j);
    end
end

p_joint_matrix = p_joint_matrix/sum(sum(p_joint_matrix));

p_test = sum(sum(triu(p_joint_matrix)));

end

