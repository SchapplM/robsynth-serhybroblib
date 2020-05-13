% Jacobian of implicit kinematic constraints of
% fourbarprisIC
% with respect to passive joint coordinates
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
% Phi_p [(no of constraints)x(no of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_p = fourbarprisIC_kinconstr_impl_pas_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:30
% EndTime: 2020-05-07 09:59:30
% DurationCPUTime: 0.01s
% Computational Cost: add. (3->3), mult. (4->4), div. (0->0), fcn. (4->4), ass. (0->2)
t4 = pkin(3) + qJ(2);
t1 = [sin(qJ(1)) * t4, -sin(qJ(3)) * pkin(2), 0; -t4 * cos(qJ(1)), cos(qJ(3)) * pkin(2), 0; -1, 1, 1;];
Phi_p = t1;
