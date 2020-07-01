% Calculate minimal parameter regressor of joint inertia matrix for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = palh2m2DE_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:46
% EndTime: 2020-06-30 17:56:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (27->7), mult. (61->32), div. (0->0), fcn. (55->6), ass. (0->21)
t13 = cos(qJ(3));
t22 = 0.2e1 * t13;
t14 = cos(qJ(2));
t21 = 0.2e1 * t14;
t16 = t14 * pkin(4) + pkin(1);
t2 = pkin(2) + t16;
t15 = t13 * pkin(5) + t2;
t1 = pkin(3) + t15;
t9 = sin(qJ(4));
t20 = t9 * t1;
t11 = sin(qJ(2));
t19 = pkin(4) * t11;
t10 = sin(qJ(3));
t18 = pkin(5) * t10;
t12 = cos(qJ(4));
t17 = t12 * t1;
t6 = t12 * t18;
t5 = t12 * t19;
t4 = t9 * t19;
t3 = t9 * t18;
t7 = [1, 0, 0, t11 ^ 2, t11 * t21, 0, 0, 0, pkin(1) * t21, -0.2e1 * pkin(1) * t11, 0.2e1 * t16, t10 ^ 2, t10 * t22, 0, 0, 0, t2 * t22, -0.2e1 * t10 * t2, 0.2e1 * t15, 1, 0.2e1 * t17, -0.2e1 * t20; 0, 0, 0, 0, 0, t11, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t13, 0, 0, 0, 0, 0, t3, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pkin(4) * (t11 * t10 + t14 * t13), pkin(4) * (-t14 * t10 + t13 * t11), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t17, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t7;
