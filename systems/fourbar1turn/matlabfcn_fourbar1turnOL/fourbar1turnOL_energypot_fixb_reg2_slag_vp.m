% Calculate inertial parameters regressor of potential energy for
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
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar1turnOL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:41:04
% EndTime: 2020-04-12 19:41:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (31->18), mult. (55->19), div. (0->0), fcn. (47->8), ass. (0->14)
t35 = g(3) * pkin(5);
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t34 = g(1) * t32 + g(2) * t29;
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t33 = -g(3) * t28 - t34 * t31;
t30 = cos(qJ(4));
t27 = sin(qJ(4));
t26 = qJ(2) + qJ(3);
t24 = cos(t26);
t23 = sin(t26);
t22 = g(1) * t29 - g(2) * t32;
t1 = [0, 0, 0, 0, 0, 0, -t34, t22, -g(3), -t35, 0, 0, 0, 0, 0, 0, t33, -g(3) * t31 + t34 * t28, -t22, -t35, 0, 0, 0, 0, 0, 0, g(3) * t23 + t34 * t24, g(3) * t24 - t34 * t23, -t22, t33 * pkin(2) - t35, 0, 0, 0, 0, 0, 0, -g(3) * t27 - t34 * t30, -g(3) * t30 + t34 * t27, -t22, -t34 * pkin(1) - t35;];
U_reg = t1;
