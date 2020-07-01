% Calculate minimal parameter regressor of potential energy for
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
% U_reg [1x9]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:30
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbarprisOL_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_energypot_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:29:57
% EndTime: 2020-06-27 17:29:57
% DurationCPUTime: 0.02s
% Computational Cost: add. (8->7), mult. (16->11), div. (0->0), fcn. (14->4), ass. (0->8)
t6 = sin(qJ(1));
t9 = g(2) * t6;
t8 = cos(qJ(1));
t7 = cos(qJ(3));
t5 = sin(qJ(3));
t4 = g(1) * t8 + t9;
t3 = g(1) * t6 - g(2) * t8;
t1 = [0, t4, -t3, t3, t4, -g(1) * (-t8 * qJ(2) + pkin(1)) + qJ(2) * t9, 0, g(1) * t7 + g(2) * t5, -g(1) * t5 + g(2) * t7;];
U_reg = t1;
