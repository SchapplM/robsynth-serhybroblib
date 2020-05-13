% Explicit kinematic constraints of
% fourbar2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% 
% Output:
% jv [5x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:17
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = fourbar2TE_kinconstr_expl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2TE_kinconstr_expl_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2TE_kinconstr_expl_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:17:20
% EndTime: 2020-04-24 20:17:20
% DurationCPUTime: 0.02s
% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (9->4), ass. (0->4)
t3 = cos(qJ(1));
t2 = sin(qJ(1));
t1 = atan2(t2, t3);
t4 = [qJ(1); atan2(t2, -t3); t1; t1; 0;];
jv = t4(:);
