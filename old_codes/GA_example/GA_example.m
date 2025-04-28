% File: ga_example.m
% Genetic Algorithm to maximize f(x) = x * sin(10 * pi * x) + 2 for x in [0, 1].

clc;
clear;

% Parameters
pop_size = 20;       % Population size
max_gen = 1000;        % Maximum number of generations
mutation_rate = 0.01; % Mutation probability
crossover_rate = 0.8; % Crossover probability
var_bounds = [0, 1]; % Variable bounds [lower, upper]

% Fitness function
fitness_func = @(x) x .* sin(10 * pi * x) + 2.0;

% Initialize population
population = var_bounds(1) + (var_bounds(2) - var_bounds(1)) * rand(pop_size, 1);

% Evolution loop
for gen = 1:max_gen
    % Evaluate fitness
    fitness = fitness_func(population);
    
    % Selection: Roulette Wheel
    fitness_sum = sum(fitness);
    prob = fitness / fitness_sum;
    cum_prob = cumsum(prob);
    new_population = zeros(pop_size, 1);
    for i = 1:pop_size
        r = rand;
        idx = find(cum_prob >= r, 1, 'first');
        new_population(i) = population(idx);
    end
    population = new_population;
    
    % Crossover
    for i = 1:2:pop_size
        if rand < crossover_rate
            p1 = population(i);
            p2 = population(i+1);
            alpha = rand; % Crossover factor
            offspring1 = alpha * p1 + (1 - alpha) * p2;
            offspring2 = alpha * p2 + (1 - alpha) * p1;
            population(i) = offspring1;
            population(i+1) = offspring2;
        end
    end
    
    % Mutation
    for i = 1:pop_size
        if rand < mutation_rate
            mutation = (var_bounds(2) - var_bounds(1)) * (rand - 0.5);
            population(i) = population(i) + mutation;
            % Ensure within bounds
            population(i) = max(var_bounds(1), min(var_bounds(2), population(i)));
        end
    end
    
    % Display best fitness in current generation
    best_fitness = max(fitness_func(population));
    fprintf('Generation %d: Best Fitness = %.4f\n', gen, best_fitness);
end

% Display final result
final_population = population;
final_fitness = fitness_func(final_population);
[~, idx] = max(final_fitness);
best_solution = final_population(idx);
fprintf('\nOptimal Solution: x = %.4f, f(x) = %.4f\n', best_solution, max(final_fitness));
