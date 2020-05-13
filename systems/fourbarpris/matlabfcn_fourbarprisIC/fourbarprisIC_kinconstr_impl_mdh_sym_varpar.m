% Implicit kinematic constraints of
% fourbarprisIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% h [(no of constraints)x1]
%   Implicit constraint equations (e.g. closed loops)
%   In a valid robot configuration qJ this has to be zero

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function h = fourbarprisIC_kinconstr_impl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_kinconstr_impl_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_kinconstr_impl_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:30
% EndTime: 2020-05-07 09:59:30
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->7), mult. (6->4), div. (0->0), fcn. (6->4), ass. (0->2)
t3 = -pkin(3) - qJ(2);
t1 = [cos(qJ(3)) * pkin(2) + pkin(1) + t3 * cos(qJ(1)); sin(qJ(3)) * pkin(2) + t3 * sin(qJ(1)); pi - qJ(1) + qJ(3) + qJ(4);];
h = t1;
