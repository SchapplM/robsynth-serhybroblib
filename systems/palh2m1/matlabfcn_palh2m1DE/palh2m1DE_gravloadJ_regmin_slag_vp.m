% Calculate minimal parameter regressor of gravitation load for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% taug_reg [4x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh2m1DE_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:35
% EndTime: 2020-06-30 17:39:36
% DurationCPUTime: 0.09s
% Computational Cost: add. (46->14), mult. (96->28), div. (0->0), fcn. (98->8), ass. (0->17)
t17 = cos(qJ(1));
t13 = sin(qJ(1));
t8 = g(1) * t17 + g(2) * t13;
t16 = cos(qJ(2));
t15 = cos(qJ(3));
t14 = cos(qJ(4));
t12 = sin(qJ(2));
t11 = sin(qJ(3));
t10 = sin(qJ(4));
t7 = g(1) * t13 - g(2) * t17;
t6 = -g(3) * t11 + t15 * t8;
t5 = g(3) * t15 + t11 * t8;
t4 = t10 * t8 + t14 * t7;
t3 = -t10 * t7 + t14 * t8;
t2 = t12 * t6 + t16 * t5;
t1 = -t12 * t5 + t16 * t6;
t9 = [0, t7, t8, 0, 0, 0, 0, 0, t7 * t16, -t12 * t7, 0, 0, 0, 0, 0, t7 * (-t11 * t12 + t15 * t16), -t7 * (t11 * t16 + t12 * t15), t7, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t16 + t12 * t8, -g(3) * t12 + t16 * t8, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3;];
taug_reg = t9;
