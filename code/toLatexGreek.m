function latexStr = toLatexGreek(name)
% toLatexGreek  Wandelt z.B. 'beta1' in '$\beta_1$' um
%
% Beispiel:
%   toLatexGreek('beta1')   ->  '$\beta_1$'
%   toLatexGreek('alpha12') ->  '$\alpha_{12}$'
%   toLatexGreek('gamma')   ->  '$\gamma$'

    % Sicherstellen, dass der Input ein String ist
    name = char(name);

    % Liste griechischer Buchstaben (klein)
    greek = ["alpha","beta","gamma","delta","epsilon","zeta","eta","theta", ...
             "iota","kappa","lambda","mu","nu","xi","omicron","pi","rho", ...
             "sigma","tau","upsilon","phi","chi","psi","omega"];

    latexStr = name; % Default: unverändert
    
    % Prüfen, ob der Name mit einem griechischen Buchstaben beginnt
    for g = greek
        if startsWith(name, g, 'IgnoreCase', true)
            rest = extractAfter(name, strlength(g)); % alles nach dem Buchstaben
            if isempty(rest)
                latexStr = ['$\' char(g) '$'];           % nur der Buchstabe
            elseif all(isstrprop(rest,'digit'))
                if strlength(rest) == 1
                    latexStr = ['$\' char(g) '_' char(rest) '$']; % z.B. β₁
                else
                    latexStr = ['$\' char(g) '_{' char(rest) '}$']; % z.B. β₁₂
                end
            else
                latexStr = ['$\' char(g) char(rest) '$'];      % z.B. βhat o.ä.
            end
            return;
        end
    end

    % Falls kein griechischer Name erkannt wurde
    latexStr = name;
end
