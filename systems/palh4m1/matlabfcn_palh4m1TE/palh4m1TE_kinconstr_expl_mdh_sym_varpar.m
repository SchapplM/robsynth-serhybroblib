% Explicit kinematic constraints of
% palh4m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AD,CB,CE,EP,HC,OT,TA,TD]';
% 
% Output:
% jv [9x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 21:48
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = palh4m1TE_kinconstr_expl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1TE_kinconstr_expl_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1TE_kinconstr_expl_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-11 21:35:19
% EndTime: 2020-04-11 21:35:19
% DurationCPUTime: 0.30s
% Computational Cost: add. (131->24), mult. (133->35), div. (12->3), fcn. (40->7), ass. (0->28)
t16 = pkin(1) ^ 2;
t22 = pkin(2) ^ 2 + t16;
t11 = sin(qJ(5));
t26 = t11 * pkin(1);
t21 = pkin(2) * t26;
t7 = -0.2e1 * t21;
t18 = t7 + t22;
t4 = 0.1e1 / t18;
t10 = qJ(2) + pkin(6);
t9 = 0.1e1 / t10;
t28 = t4 * t9;
t19 = pkin(3) - t10;
t20 = -pkin(3) - t10;
t23 = t16 + t7;
t1 = sqrt(-((pkin(2) - t20) * (pkin(2) + t20) + t23) * ((pkin(2) - t19) * (pkin(2) + t19) + t23));
t12 = cos(qJ(5));
t27 = t1 * t12;
t14 = 0.1e1 / pkin(3);
t25 = t14 * t4;
t24 = t14 * t9;
t13 = pkin(3) ^ 2;
t17 = -t13 + t22;
t8 = t10 ^ 2;
t6 = -pkin(2) * t11 + pkin(1);
t5 = -pkin(2) + t26;
t3 = t13 - t8 + t18;
t2 = t7 + t8 + t17;
t15 = [qJ(1); atan2((pkin(1) * t12 * t2 - t5 * t1) * t28, (-pkin(1) * t27 - t5 * t2) * t28); t10; atan2(t1 * t24, (-t17 + t8 + 0.2e1 * t21) * t24); qJ(3); qJ(4); qJ(5); atan2((pkin(2) * t12 * t3 + t6 * t1) * t25 / 0.2e1, -(-pkin(2) * t27 + t6 * t3) * t25 / 0.2e1); 0;];
jv = t15(:);
