% Calculate minimal parameter regressor of potential energy for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:56
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar1turnOL_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:56:28
% EndTime: 2020-06-27 16:56:28
% DurationCPUTime: 0.03s
% Computational Cost: add. (20->10), mult. (34->16), div. (0->0), fcn. (34->8), ass. (0->11)
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t29 = g(1) * t28 + g(2) * t25;
t27 = cos(qJ(2));
t26 = cos(qJ(4));
t24 = sin(qJ(2));
t23 = sin(qJ(4));
t22 = qJ(2) + qJ(3);
t21 = cos(t22);
t20 = sin(t22);
t1 = [0, -t29, g(1) * t25 - g(2) * t28, 0, 0, 0, 0, 0, -g(3) * t24 - t29 * t27, -g(3) * t27 + t29 * t24, 0, 0, 0, 0, 0, g(3) * t20 + t29 * t21, g(3) * t21 - t29 * t20, 0, 0, 0, 0, 0, -g(3) * t23 - t29 * t26, -g(3) * t26 + t29 * t23;];
U_reg = t1;
