% Calculate inertial parameters regressor of potential energy for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbarprisOL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_energypot_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:07
% EndTime: 2020-05-07 09:52:07
% DurationCPUTime: 0.03s
% Computational Cost: add. (11->10), mult. (17->12), div. (0->0), fcn. (14->4), ass. (0->8)
t6 = sin(qJ(1));
t9 = g(2) * t6;
t8 = cos(qJ(1));
t7 = cos(qJ(3));
t5 = sin(qJ(3));
t4 = g(1) * t8 + t9;
t3 = g(1) * t6 - g(2) * t8;
t1 = [0, 0, 0, 0, 0, 0, t4, -t3, -g(3), -g(1) * pkin(1), 0, 0, 0, 0, 0, 0, t3, g(3), t4, -g(1) * (-t8 * qJ(2) + pkin(1)) + qJ(2) * t9, 0, 0, 0, 0, 0, 0, g(1) * t7 + g(2) * t5, -g(1) * t5 + g(2) * t7, -g(3), 0;];
U_reg = t1;
