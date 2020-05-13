% Explicit kinematic constraints of
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% jv [6x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = fourbar1turnTE_kinconstr_expl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_kinconstr_expl_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_kinconstr_expl_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:23
% EndTime: 2020-04-12 19:18:23
% DurationCPUTime: 0.24s
% Computational Cost: add. (101->23), mult. (133->35), div. (12->3), fcn. (40->7), ass. (0->27)
t9 = cos(qJ(2));
t27 = pkin(2) * t9;
t15 = pkin(1) ^ 2;
t18 = pkin(1) * t27;
t7 = -0.2e1 * t18;
t21 = t15 + t7;
t24 = -pkin(3) + pkin(4);
t25 = -pkin(3) - pkin(4);
t1 = sqrt(-((pkin(2) - t25) * (pkin(2) + t25) + t21) * ((pkin(2) - t24) * (pkin(2) + t24) + t21));
t8 = sin(qJ(2));
t26 = t1 * t8;
t11 = 0.1e1 / pkin(4);
t19 = pkin(2) ^ 2 + t15;
t17 = t7 + t19;
t4 = 0.1e1 / t17;
t23 = t11 * t4;
t13 = 0.1e1 / pkin(3);
t22 = t13 * t4;
t20 = t11 * t13;
t12 = pkin(3) ^ 2;
t16 = -t12 + t19;
t10 = pkin(4) ^ 2;
t6 = pkin(1) * t9 - pkin(2);
t5 = pkin(1) - t27;
t3 = t10 + t7 + t16;
t2 = -t10 + t12 + t17;
t14 = [qJ(1); qJ(2); atan2((pkin(1) * t8 * t2 - t6 * t1) * t22, (-pkin(1) * t26 - t6 * t2) * t22); atan2((pkin(2) * t8 * t3 + t5 * t1) * t23 / 0.2e1, -(-pkin(2) * t26 + t5 * t3) * t23 / 0.2e1); atan2(t1 * t20, (t10 - t16 + 0.2e1 * t18) * t20); 0;];
jv = t14(:);
