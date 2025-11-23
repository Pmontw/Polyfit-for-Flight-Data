clf;clc;clear
load count.dat
data = readtable("Full Test Flight negative.csv");
data_start = 1;
data_end = 24000;
first_column = 10; 

% Change next 6 values for desired power of plynomial
rotx_power = 6;
roty_power = 6;
rotz_power = 6;

accx_power = 37;
accy_power = 35;
accz_power = 37;

% Array for each component of rotation and acceleration in .01s table
rot(data_start:data_end,1) = table2array(data(data_start:data_end,first_column));
rot(data_start:data_end,2) = table2array(data(data_start:data_end,first_column+1));
rot(data_start:data_end,3) = table2array(data(data_start:data_end,first_column+2));
acc(data_start:data_end,1) = table2array(data(data_start:data_end,first_column+3));
acc(data_start:data_end,2) = table2array(data(data_start:data_end,first_column+4));
acc(data_start:data_end,3) = table2array(data(data_start:data_end,first_column+5));

% varaible for each array of data
rotx_data= rot(data_start:data_end,1);
roty_data= rot(data_start:data_end,2);
rotz_data= rot(data_start:data_end,3);
accx_data= acc(data_start:data_end,1);
accy_data= acc(data_start:data_end,2);
accz_data= acc(data_start:data_end,3);

%Time for each graph of the components
sample_rate = 0.001;  % seconds per sample
t = ((data_start:data_end) - 1) * sample_rate;

Time1=t;
Time2=t;
Time3=t;
Time4=t;
Time5=t;
Time6=t;

direction_array = [1,2,3]; %array to represent x,y,z



%Polynimial functions of the rotational data
function_rot_array = [];

% Fit and evaluate rotx
[rotx, ~, mu_rotx] = polyfit(Time1, rot(data_start:data_end,1), rotx_power);
rotx2 = polyval(rotx, (Time1 - mu_rotx(1)) / mu_rotx(2));
function_rot_array = [function_rot_array, rotx2(:)];

% Fit and evaluate roty
[roty, ~, mu_roty] = polyfit(Time2, rot(data_start:data_end,2), roty_power);
roty2 = polyval(roty, (Time2 - mu_roty(1)) / mu_roty(2));
function_rot_array = [function_rot_array, roty2(:)];

% Fit and evaluate rotz
[rotz, ~, mu_rotz] = polyfit(Time3, rot(data_start:data_end,3), rotz_power);
rotz2 = polyval(rotz, (Time3 - mu_rotz(1)) / mu_rotz(2));
function_rot_array = [function_rot_array, rotz2(:)];



%plot for rotations 
rot_names = ["Rot X", "Rot Y", "Rot Z"];
rot_equations = {rotx, roty, rotz};
rot_mu = {mu_rotx, mu_roty, mu_rotz};

for i = direction_array
    figure
    hold on
    plot(t, rot(data_start:data_end, i), 'o', 'DisplayName', 'Raw Data');
    plot(t, function_rot_array(:, i), '-', 'LineWidth', 2, 'DisplayName', 'Polynomial Fit');
    grid on;
    
    % Compute R²
    y_actual = rot(data_start:data_end,i);
    y_predicted = function_rot_array(:,i);
    R2 = 1 - sum((y_actual - y_predicted).^2) / sum((y_actual - mean(y_actual)).^2);

    % Generate symbolic equation
    syms t_sym x
    coeffs = rot_equations{i};
    mu_tmp = rot_mu{i};
    Eq_sym = poly2sym(coeffs, x);
    Eq_unscaled = expand(subs(Eq_sym, x, (t_sym - mu_tmp(1)) / mu_tmp(2)));
    eq_str = char(vpa(Eq_unscaled, 3)); % shorten precision for readability

    % Title with R² and equation
    title({rot_names(i), ...
        sprintf('R² = %.4f', R2)});

    legend('Location', 'best');
    xlabel('Time [s]');
    ylabel('[rad/s]');
    hold off
end






function_acc_array = [];

[accx, ~, mu_accx] = polyfit(Time4, acc(data_start:data_end,1), accx_power); % Fit with centering & scaling

accx2 = polyval(accx, (Time4 - mu_accx(1)) / mu_accx(2)); % Apply same transformation to accx
function_acc_array = [function_acc_array, accx2(:)];

[accy, ~, mu] = polyfit(Time5, acc(data_start:data_end,2), accy_power); % Fit with centering & scaling

accy2 = polyval(accy, (Time5 - mu(1)) / mu(2)); % Apply same transformation to accy
function_acc_array = [function_acc_array, accy2(:)];

[accz, ~, mu] = polyfit(Time6, acc(data_start:data_end,3), accz_power); % Fit with centering & scaling

accz2 = polyval(accz, (Time6 - mu(1)) / mu(2)); % Apply same transformation to x2
function_acc_array = [function_acc_array, accz2(:)];

%plot for acclerations 
acc_names = ["Acc X", "Acc Y", "Acc Z"];
acc_equations = {accx, accy, accz};
acc_mu = {mu_accx, mu, mu}; % use appropriate mu for each if different

for i = direction_array
    figure
    hold on
    plot(t, acc(data_start:data_end, i), 'o', 'DisplayName', 'Raw Data');
    plot(t, function_acc_array(:, i), '-', 'LineWidth', 2, 'DisplayName', 'Polynomial Fit');
    grid on;

    % Compute R²
    y_actual = acc(data_start:data_end,i);
    y_predicted = function_acc_array(:,i);
    R2 = 1 - sum((y_actual - y_predicted).^2) / sum((y_actual - mean(y_actual)).^2);

    % Generate symbolic equation
    syms t_sym x
    coeffs = acc_equations{i};
    mu_tmp = acc_mu{i};
    Eq_sym = poly2sym(coeffs, x);
    Eq_unscaled = expand(subs(Eq_sym, x, (t_sym - mu_tmp(1)) / mu_tmp(2)));
    eq_str = char(vpa(Eq_unscaled, 3)); % shorten precision for readability

    % Title with R² and equation
    title({acc_names(i), ...
        sprintf('R² = %.4f', R2)});
        

    legend('Location', 'best');
    xlabel('Time [s]');
    ylabel('[m/s²]');
    hold off
end



% Compute R² of rotation

for direction = direction_array
    y_actual = rot(data_start:data_end,direction);  % Original acceleration data
    y_predicted = function_rot_array(:,direction);       % Polynomial fitted values

    y_actual = y_actual(:);  % Force as column vector if necessary
    y_predicted = y_predicted(:);  % Force as column vector if necessary
    
    residuals = y_actual - y_predicted;
    
    SS_res = sum(residuals.^2);  
    SS_tot = sum((y_actual - mean(y_actual)).^2); 
    
    R2 = 1 - (SS_res/SS_tot);
    
    % Display R² in command window
    fprintf('Rot R² value: %.4f\n', R2);

end

%display rotation equations
rot_name = ["Rotx","Roty","Rotz"];

% Store polynomial coefficients in a struct
rot_coeffs.Rotx = rotx;
rot_coeffs.Roty = roty;
rot_coeffs.Rotz = rotz;

% Iterate through each rotation name
for i = 1:length(rot_name)
    syms t x
    name = rot_name(i); % Get string name
    coeffs = rot_coeffs.(char(name)); % Access the corresponding coefficients
   
    % Convert polynomial to symbolic equation
    EqRot_norm = poly2sym(coeffs, x); % Step 1: Create symbolic polynomial in x (normalized variable)
    norm = (t - mu(1)) / mu(2);
   
    EqRot_unscaled = expand(subs(EqRot_norm, x, norm));
    
    % Display results
    disp(name);  % Display the name (e.g., Rotx")
    
    % Extract the coefficients and powers of t from the symbolic equation
    terms = sym2poly(EqRot_unscaled); % Extract terms of the equation
    powers = fliplr(0:length(terms)-1); % Get the powers of t for the polynomial
    
    % Loop through the terms and format the output
    output_str = '';
    UDF_Output_Stur =  ''; 
    for j = 1:length(terms)
        coeff = terms(j);
        power = powers(j);
        
        % Format the coefficient in scientific notation with 4 decimal places
        coeff_str = sprintf('%1.4e', coeff);
        
        % Break coefficient into base and exponent
        [coeff_base, coeff_exp] = strtok(coeff_str, 'e');
        coeff_exp = str2double(coeff_exp(2:end)); % Get exponent
        
        % Build the formatted string for this term
        term_str = sprintf('%1.4f [rad/s^%d] * 10^%d * t^%d', str2double(coeff_base), power + 1, coeff_exp, power);
        UDF_str = sprintf('%1.4f * pow(10.0,%d.0) * pow(time,%d)', str2double(coeff_base), coeff_exp, power); 
        
        
        % Append to the output string
        output_str = strcat(output_str, term_str);
        UDF_Output_Stur = strcat(UDF_Output_Stur, UDF_str);
        % Add a "+" sign except for the last term
        if j < length(terms)
            output_str = strcat(output_str, ' + ');
            UDF_Output_Stur = strcat(UDF_Output_Stur, ' + ');
        end
    end
    
    % Display the output string in the desired format
    disp(output_str);
    disp(UDF_Output_Stur)
end


% Compute R² of acceleration

for direction = direction_array
    y_actual = acc(data_start:data_end,direction);  % Original acceleration data
    y_predicted = function_acc_array(:,direction);       % Polynomial fitted values

    y_actual = y_actual(:);  % Force as column vector if necessary
    y_predicted = y_predicted(:);  % Force as column vector if necessary
    
    residuals = y_actual - y_predicted;
    
    SS_res = sum(residuals.^2);  
    SS_tot = sum((y_actual - mean(y_actual)).^2); 
    
    R2 = 1 - (SS_res/SS_tot);
    
    % Display R² in command window
    fprintf('Acc R² value: %.4f\n', R2);

end



%display of y acceleration equation
acc_name = ["Accx","Accy","Accz"];

% Store polynomial coefficients in a struct
acc_coeffs.Accx = accx;
acc_coeffs.Accy = accy;
acc_coeffs.Accz = accz;


% Iterate through each rotation name
for i = 1:length(acc_name)
    syms t x
    name = acc_name(i); % Get string name
    coeffs = acc_coeffs.(char(name)); % Access the corresponding coefficients
   
    % Convert polynomial to symbolic equation
    EqAcc_norm = poly2sym(coeffs, x); % Step 1: Create symbolic polynomial in x (normalized variable)
    norm = (t - mu(1)) / mu(2);
   
    EqAcc_unscaled = expand(subs(EqAcc_norm, x, norm));
    
    % Display results
    disp(name);  % Display the name (e.g., Rotx")
    
    % Extract the coefficients and powers of t from the symbolic equation
    terms = sym2poly(EqAcc_unscaled); % Extract terms of the equation
    powers = fliplr(0:length(terms)-1); % Get the powers of t for the polynomial
    
    % Loop through the terms and format the output
    output_str = '';
    UDF_Output_Stur =  ''; 
    for j = 1:length(terms)
        coeff = terms(j);
        power = powers(j);
        
        % Format the coefficient in scientific notation with 4 decimal places
        coeff_str = sprintf('%1.4e', coeff);
        
        % Break coefficient into base and exponent
        [coeff_base, coeff_exp] = strtok(coeff_str, 'e');
        coeff_exp = str2double(coeff_exp(2:end)); % Get exponent
        
        % Build the formatted string for this term
        term_str = sprintf('%1.4f [m/s^%d] * 10^%d * t^%d', str2double(coeff_base), power + 2, coeff_exp, power);
        UDF_str = sprintf('%1.4f * pow(10,%d) * pow(time,%d)', str2double(coeff_base), coeff_exp, power); 
        % Append to the output string
        output_str = strcat(output_str, term_str);
        UDF_Output_Stur = strcat(UDF_Output_Stur, UDF_str);
        % Add a "+" sign except for the last term
        if j < length(terms)
            output_str = strcat(output_str, ' + ');
            UDF_Output_Stur = strcat(UDF_Output_Stur, ' + ');
        end
    end
    
    % Display the output string in the desired format
    disp(output_str);
    disp(UDF_Output_Stur)
end

function display_poly_equation(time, data, degree, var_name, unit_power)
    % Fit polynomial
    [p, ~, mu] = polyfit(time, data, degree);

    % Symbolic expression
    syms x t
    poly_norm = poly2sym(p, x);

    % Substitute and expand
    x_to_t = (t - mu(1)) / mu(2);
    poly_unscaled = expand(subs(poly_norm, x, x_to_t));

    % Extract coefficients and powers
    terms = sym2poly(poly_unscaled);
    powers = fliplr(0:length(terms)-1);

    % Display normalized form
    fprintf('\n%s (normalized):\n', var_name);
    fprintf('f(x) = %s\n', char(vpa(poly_norm, 6)));
    fprintf('where x = (t - %.6f) / %.6f\n\n', mu(1), mu(2));

    % Display unscaled form
    fprintf('%s (expanded in time):\n', var_name);
    disp(vpa(poly_unscaled, 8));

    % Construct formatted string outputs
    output_str = '';
    UDF_Output_Str = '';
    for j = 1:length(terms)
        coeff = terms(j);
        power = powers(j);

        coeff_str = sprintf('%1.4e', coeff);
        [coeff_base, coeff_exp] = strtok(coeff_str, 'e');
        coeff_exp = str2double(coeff_exp(2:end));

        term_str = sprintf('%1.4f [%s^%d] * 10^%d * t^%d', str2double(coeff_base), unit_power.unit, power + unit_power.offset, coeff_exp, power);
        UDF_str = sprintf('%1.4f * pow(10.0,%d.0) * pow(time,%d)', str2double(coeff_base), coeff_exp, power);

        output_str = strcat(output_str, term_str);
        UDF_Output_Str = strcat(UDF_Output_Str, UDF_str);

        if j < length(terms)
            output_str = strcat(output_str, ' + ');
            UDF_Output_Str = strcat(UDF_Output_Str, ' + ');
        end
    end

    % Display final formatted strings
    fprintf('%s formatted:\n%s\n', var_name, output_str);
    fprintf('%s UDF format:\n%s\n', var_name, UDF_Output_Str);
end

