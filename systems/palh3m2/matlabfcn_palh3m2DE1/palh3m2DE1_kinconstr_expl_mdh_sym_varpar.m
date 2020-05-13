% Explicit kinematic constraints of
% palh3m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% jv [12x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = palh3m2DE1_kinconstr_expl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_kinconstr_expl_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_kinconstr_expl_mdh_sym_varpar: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:57:18
% EndTime: 2020-05-07 01:57:18
% DurationCPUTime: 0.04s
% Computational Cost: add. (88->14), mult. (144->24), div. (0->0), fcn. (234->19), ass. (0->23)
t23 = cos(pkin(15));
t22 = cos(qJ(2));
t21 = cos(pkin(14));
t20 = cos(qJ(3));
t19 = sin(pkin(14));
t18 = sin(pkin(15));
t17 = sin(qJ(2));
t16 = sin(qJ(3));
t15 = cos(pkin(18));
t14 = sin(pkin(18));
t13 = cos(pkin(16));
t12 = sin(pkin(16));
t11 = pkin(17) + pkin(18);
t10 = cos(t11);
t9 = sin(t11);
t7 = t17 * t23 + t22 * t18;
t6 = t17 * t18 - t22 * t23;
t5 = -t6 * t16 + t7 * t20;
t4 = t7 * t16 + t6 * t20;
t3 = t10 * t4 + t9 * t5;
t2 = -t10 * t5 + t9 * t4;
t1 = atan2(t2, t3);
t8 = [qJ(1); qJ(2); qJ(3); atan2(t3 * t12 + t2 * t13, t2 * t12 - t3 * t13); qJ(4); atan2(t6 * t19 + t7 * t21, t7 * t19 - t6 * t21); atan2(-t6 * t14 + t7 * t15, t7 * t14 + t6 * t15); t1; atan2(t7, -t6); t1; 0; 0;];
jv = t8(:);
