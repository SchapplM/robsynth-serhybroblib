% Calculate minimal parameter regressor of gravitation load for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% taug_reg [4x9]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:30
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbarprisOL_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_gravloadJ_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:29:58
% EndTime: 2020-06-27 17:29:58
% DurationCPUTime: 0.06s
% Computational Cost: add. (8->5), mult. (18->9), div. (0->0), fcn. (16->4), ass. (0->7)
t4 = sin(qJ(1));
t6 = cos(qJ(1));
t1 = -g(1) * t4 + g(2) * t6;
t5 = cos(qJ(3));
t3 = sin(qJ(3));
t2 = g(1) * t6 + g(2) * t4;
t7 = [0, t1, -t2, t2, t1, t1 * qJ(2), 0, 0, 0; 0, 0, 0, 0, 0, t2, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t5, -g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t7;
