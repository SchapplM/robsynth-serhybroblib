% Calculate minimal parameter regressor of gravitation load for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% taug_reg [4x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh2m2DE_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:46
% EndTime: 2020-06-30 17:56:46
% DurationCPUTime: 0.07s
% Computational Cost: add. (28->10), mult. (64->20), div. (0->0), fcn. (64->8), ass. (0->13)
t13 = cos(qJ(1));
t9 = sin(qJ(1));
t4 = g(1) * t13 + g(2) * t9;
t12 = cos(qJ(2));
t11 = cos(qJ(3));
t10 = cos(qJ(4));
t8 = sin(qJ(2));
t7 = sin(qJ(3));
t6 = sin(qJ(4));
t3 = g(1) * t9 - g(2) * t13;
t2 = t3 * t10 + t6 * t4;
t1 = t10 * t4 - t6 * t3;
t5 = [0, t3, t4, 0, 0, 0, 0, 0, t3 * t12, -t8 * t3, t3, 0, 0, 0, 0, 0, t3 * t11, -t7 * t3, t3, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, -t12 * g(3) + t4 * t8, g(3) * t8 + t4 * t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * g(3) + t4 * t7, g(3) * t7 + t4 * t11, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t5;
