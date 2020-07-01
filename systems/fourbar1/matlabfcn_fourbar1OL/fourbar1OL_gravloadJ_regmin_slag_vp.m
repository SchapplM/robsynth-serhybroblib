% Calculate minimal parameter regressor of gravitation load for
% fourbar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% taug_reg [4x9]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:43
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar1OL_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1OL_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_gravloadJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:43:19
% EndTime: 2020-06-26 17:43:19
% DurationCPUTime: 0.06s
% Computational Cost: add. (16->7), mult. (16->12), div. (0->0), fcn. (16->6), ass. (0->10)
t9 = cos(qJ(1));
t8 = cos(qJ(3));
t7 = sin(qJ(1));
t6 = sin(qJ(3));
t5 = qJ(1) + qJ(2);
t4 = cos(t5);
t3 = sin(t5);
t2 = -g(1) * t4 - g(2) * t3;
t1 = -g(1) * t3 + g(2) * t4;
t10 = [0, g(1) * t7 - g(2) * t9, g(1) * t9 + g(2) * t7, 0, t1, t2, 0, 0, 0; 0, 0, 0, 0, t1, t2, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, g(1) * t6 - g(2) * t8, g(1) * t8 + g(2) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t10;
